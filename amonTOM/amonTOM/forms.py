from django import forms
from astropy import units as u
from django.forms import ValidationError
from django.conf import settings
from django.contrib.auth.models import Group
from guardian.shortcuts import assign_perm, get_groups_with_perms, remove_perm
from tom_observations.models import ObservationRecord

from tom_targets.models import (
    Target, TargetExtra, SIDEREAL_FIELDS, NON_SIDEREAL_FIELDS, REQUIRED_SIDEREAL_FIELDS,
    REQUIRED_NON_SIDEREAL_FIELDS, REQUIRED_NON_SIDEREAL_FIELDS_PER_SCHEME
)
from tom_targets.forms import CoordinateField, extra_field_to_form_field
import logging

logger = logging.getLogger(__name__)

class TargetForm(forms.ModelForm):
    groups = forms.ModelMultipleChoiceField(Group.objects.none(), required=True, widget=forms.CheckboxSelectMultiple)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extra_fields = {}
        for extra_field in settings.EXTRA_FIELDS:
            # Add extra fields to the form
            field_name = extra_field['name']
            self.extra_fields[field_name] = extra_field_to_form_field(extra_field['type'])
            # Populate them with initial values if this is an update
            if kwargs['instance']:
                te = TargetExtra.objects.filter(target=kwargs['instance'], key=field_name)
                if te.exists():
                    self.extra_fields[field_name].initial = te.first().typed_value(extra_field['type'])

            self.fields.update(self.extra_fields)

    def save(self, commit=True):
        instance = super().save(commit=commit)
        if commit:
            for field in settings.EXTRA_FIELDS:
                if self.cleaned_data.get(field['name']) is not None:
                    TargetExtra.objects.update_or_create(
                            target=instance,
                            key=field['name'],
                            defaults={'value': self.cleaned_data[field['name']]}
                    )
            if not 'groups' in self.cleaned_data:
                self.cleaned_data['groups'] = get_groups_with_perms(instance)
            # Save groups for this target
            for group in self.cleaned_data['groups']:
                assign_perm('tom_targets.view_target', group, instance)
                assign_perm('tom_targets.change_target', group, instance)
                assign_perm('tom_targets.delete_target', group, instance)
            for group in get_groups_with_perms(instance):
                if group not in self.cleaned_data['groups']:
                    remove_perm('tom_targets.view_target', group, instance)
                    remove_perm('tom_targets.change_target', group, instance)
                    remove_perm('tom_targets.delete_target', group, instance)

        return instance

    class Meta:
        abstract = True
        model = Target
        fields = '__all__'
        widgets = {'type': forms.HiddenInput()}

class SiderealTargetCreateForm(TargetForm):
    ra = CoordinateField(required=True, label='Right Ascension', c_type='ra',
                         help_text='Right Ascension, in decimal degrees or sexagesimal hours. See '
                                   'https://docs.astropy.org/en/stable/api/astropy.coordinates.Angle.html for '
                                   'supported sexagesimal inputs.')
    dec = CoordinateField(required=True, label='Declination', c_type='dec',
                          help_text='Declination, in decimal or sexagesimal degrees. See '
                                    ' https://docs.astropy.org/en/stable/api/astropy.coordinates.Angle.html for '
                                    'supported sexagesimal inputs.')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for field in REQUIRED_SIDEREAL_FIELDS:
            self.fields[field].required = True

    class Meta(TargetForm.Meta):
        fields = SIDEREAL_FIELDS


class NonSiderealTargetCreateForm(TargetForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for field in REQUIRED_NON_SIDEREAL_FIELDS:
            self.fields[field].required = True

    def clean(self):
        """
        Look at the 'scheme' field and check the fields required for the
        specified field have been given
        """
        cleaned_data = super().clean()
        scheme = cleaned_data['scheme']  # scheme is a required field, so this should be safe
        required_fields = REQUIRED_NON_SIDEREAL_FIELDS_PER_SCHEME[scheme]

        for field in required_fields:
            if not cleaned_data.get(field):
                # Get verbose names of required fields
                field_names = [
                    "'" + Target._meta.get_field(f).verbose_name + "'"
                    for f in required_fields
                ]
                scheme_name = dict(Target.TARGET_SCHEMES)[scheme]
                raise ValidationError(
                    "Scheme '{}' requires fields {}".format(scheme_name, ', '.join(field_names))
                )

    class Meta(TargetForm.Meta):
        fields = NON_SIDEREAL_FIELDS

class DataProductUploadForm(forms.Form):
    observation_record = forms.ModelChoiceField(
        ObservationRecord.objects.all(),
        widget=forms.HiddenInput(),
        required=False
    )
    target = forms.ModelChoiceField(
        Target.objects.all(),
        widget=forms.HiddenInput(),
        required=False
    )
    files = forms.FileField(
        widget=forms.ClearableFileInput(
            attrs={'multiple': True}
        )
    )
    data_product_type = forms.ChoiceField(
        choices=[v for k, v in settings.DATA_PRODUCT_TYPES.items()],
        widget=forms.RadioSelect(),
        required=True
    )
    referrer = forms.CharField(
        widget=forms.HiddenInput()
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not settings.TARGET_PERMISSIONS_ONLY:
            self.fields['groups'] = forms.ModelMultipleChoiceField(Group.objects.none(),
                                                                   required=True,
                                                                   widget=forms.CheckboxSelectMultiple)
