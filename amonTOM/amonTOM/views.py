import logging

import numpy as np
from django.views.generic import TemplateView, FormView, DeleteView
from django_filters.views import FilterView
from tom_alerts.views import RunQueryView, BrokerQueryListView, BrokerQueryFilter
from django.shortcuts import get_object_or_404
from tom_alerts.models import BrokerQuery
from django.utils import timezone
from django.core.cache import cache
from django.contrib.auth.mixins import LoginRequiredMixin
from guardian.mixins import PermissionListMixin
import json
from django.forms import HiddenInput
from tom_common.mixins import Raise403PermissionRequiredMixin
from tom_alerts.alerts import get_service_class, get_service_classes
from django.conf import settings
from django.urls import reverse, reverse_lazy
from django.shortcuts import redirect
from django.core.exceptions import PermissionDenied
from django.contrib.auth.forms import UserCreationForm, UsernameField
from django.contrib.auth.models import User, Group
from django import forms
from django.db import IntegrityError
from django.views.generic.detail import DetailView

from django.contrib import messages
from django.db import transaction
from django.views.generic.edit import CreateView, UpdateView
from django.views.generic import TemplateView, View

from guardian.mixins import PermissionRequiredMixin
from guardian.shortcuts import get_objects_for_user, get_groups_with_perms, assign_perm
from django.db.models import F

from tom_common.hooks import run_hook
from tom_targets.models import Target
from tom_targets.forms import TargetExtraFormset, TargetNamesFormset
from tom_observations.observing_strategy import RunStrategyForm
from .forms import SiderealTargetCreateForm, NonSiderealTargetCreateForm, DataProductUploadForm

logger = logging.getLogger(__name__)

class AboutView(TemplateView):
    template_name = 'about.html'

    def get_context_data(self, **kwargs):
        return {'targets': Target.objects.all()}

class AmonBrokerQueryListView(FilterView):
    """
    View that displays all saved ``BrokerQuery`` objects.
    """
    model = BrokerQuery
    template_name = 'tom_alerts/brokerquery_list.html'
    filterset_class = BrokerQueryFilter
    ordering = [F('last_run').desc(nulls_first=True)]
    
    def get_context_data(self, *args, **kwargs):
        """
        Adds the brokers available to the TOM to the context dictionary.
    
        :returns: context
        :rtype: dict
        """
        context = super().get_context_data(*args, **kwargs)
        private_brokers = settings.PRIVATE_BROKERS
        context['private_brokers'] = private_brokers
        context['amon_brokers'] = settings.AMON_BROKERS
        installed_brokers = get_service_classes()
        if not self.request.user.groups.filter(name='AMON members').exists():
            for broker in private_brokers:
                if broker in installed_brokers:
                    del installed_brokers[broker]

        context['installed_brokers'] = installed_brokers
        return context

class PublicAmonRunQueryView(TemplateView):
    """
    View that handles the running of a specific ``BrokerQuery``.
    """
    template_name = 'tom_alerts/query_result_amon.html'

    def get_context_data(self, *args, **kwargs):
        """
        Runs the ``fetch_alerts`` method specific to the given ``BrokerQuery`` and adds the matching alerts to the
        context dictionary.

        :returns: context
        :rtype: dict
        """
        query = get_object_or_404(BrokerQuery, pk=self.kwargs['pk'])
        if query.broker in settings.PRIVATE_BROKERS and not self.request.user.groups.filter(name='AMON members').exists():
            raise PermissionDenied()
        context = super().get_context_data()
        context['amon_brokers'] = settings.AMON_BROKERS
        context['broker'] = query.broker
        broker_class = get_service_class(query.broker)()
        alerts = broker_class.fetch_alerts(query.parameters_as_dict)
        context['alerts'] = []
        query.last_run = timezone.now()
        query.save()
        context['query'] = query
        # NB ideally I would use pagination but I didn't manage to do it
        # Anyway we won't have more than 100 public events in one broker within the coming years
        max_n_evts = 100 # maximum number of events to display
        context['max_n_evts'] = 0
        try:
            while len(context['alerts']) < max_n_evts:
                alert = next(alerts)
                generic_alert = broker_class.to_generic_alert(alert)
                cache.set('alert_{}'.format(generic_alert.id), json.dumps(alert), 3600)
                context['alerts'].append(generic_alert)
            if len(context['alerts']) == max_n_evts:
                context['max_n_evts'] = max_n_evts
        except StopIteration:
            pass
        # Reorder
        ordering = self.request.GET.get('ordering', 'defaultOrderField')
        ordering_list = []
        if ordering != 'defaultOrderField':
            for alert in context['alerts']:
                if not hasattr(alert, ordering) or type(eval('alert.'+ordering)) is type(None):
                    ordering_list += [1e10]
                else:
                    ordering_list += [eval('alert.'+ordering)]
            arg_sort = np.argsort(ordering_list)
            context['alerts'] = [context['alerts'][arg] for arg in arg_sort]
        return context

class PrivateAmonRunQueryView(LoginRequiredMixin, TemplateView):
    """
    View that handles the running of a specific ``BrokerQuery``. Requires authentication.
    """
    template_name = 'tom_alerts/query_result_amon.html'

    def get_context_data(self, *args, **kwargs):
        """
        Runs the ``fetch_alerts`` method specific to the given ``BrokerQuery`` and adds the matching alerts to the
        context dictionary.

        :returns: context
        :rtype: dict
        """
        context = super().get_context_data()
        query = get_object_or_404(BrokerQuery, pk=self.kwargs['pk'])
        context['amon_brokers'] = settings.AMON_BROKERS
        context['broker'] = query.broker
        broker_class = get_service_class(query.broker)()
        alerts = broker_class.fetch_alerts(query.parameters_as_dict)
        context['alerts'] = []
        query.last_run = timezone.now()
        query.save()
        # NB ideally I would use pagination but I didn't manage to do it
        # Anyway we won't have more than 100 public events in one broker within the coming years
        max_n_evts = 100 # maximum number of events to display
        context['max_n_evts'] = 0
        try:
            while len(context['alerts']) < max_n_evts:
                alert = next(alerts)
                generic_alert = broker_class.to_generic_alert(alert)
                cache.set('alert_{}'.format(generic_alert.id), json.dumps(alert), 3600)
                context['alerts'].append(generic_alert)
            if len(context['alerts']) == max_n_evts:
                context['max_n_evts'] = max_n_evts
        except StopIteration:
            pass
        return context


class PublicBrokerQueryDeleteView(DeleteView):
    """
    View that handles the deletion of a saved ``BrokerQuery``.
    """
    model = BrokerQuery
    success_url = reverse_lazy('tom_alerts:list')

class PublicBrokerQueryCreateView(FormView):
    """
    View for creating a new query to a broker.
    """
    template_name = 'tom_alerts/query_form.html'

    def get_broker_name(self):
        """
        Returns the broker specified in the request

        :returns: Broker name
        :rtype: str
        """
        if self.request.method == 'GET':
            broker = self.request.GET.get('broker')
        elif self.request.method == 'POST':
            broker = self.request.POST.get('broker')
        if broker in settings.PRIVATE_BROKERS and not self.request.user.groups.filter(name='AMON members').exists():
            raise PermissionDenied()
            return None
        return broker

    def get_form_class(self):
        """
        Returns the form class to use in this view. The form class will be the one defined in the specific broker
        module for which a new query is being created.
        """
        broker_name = self.get_broker_name()

        if not broker_name:
            raise ValueError('Must provide a broker name')

        return get_service_class(broker_name).form

    def get_form(self):
        """
        Returns an instance of the form to be used in this view.

        :returns: Form instance
        :rtype: django.forms.Form
        """
        form = super().get_form()
        form.helper.form_action = reverse('public-create')
        return form

    def get_initial(self):
        """
        Returns the initial data to use for forms on this view.

        :returns: dict of initial values
        :rtype: dict
        """
        initial = super().get_initial()
        initial['broker'] = self.get_broker_name()
        return initial

    def form_valid(self, form):
        """
        Saves the associated ``BrokerQuery`` and redirects to the ``BrokerQuery`` list.
        """
        form.save()
        return redirect(reverse('tom_alerts:list'))


class PublicBrokerQueryUpdateView(FormView):
    """
    View that handles the modification of a previously saved ``BrokerQuery``.
    """
    template_name = 'tom_alerts/query_form.html'

    def get_object(self):
        """
        Returns the ``BrokerQuery`` object that corresponds with the ID in the query path.

        :returns: ``BrokerQuery`` object
        :rtype: ``BrokerQuery``
        """
        query = BrokerQuery.objects.get(pk=self.kwargs['pk'])
        if query.broker in settings.PRIVATE_BROKERS and not self.request.user.groups.filter(name='AMON members').exists():
            raise PermissionDenied()
            return None
        return query

    def get_form_class(self):
        """
        Returns the form class to use in this view. The form class will be the one defined in the specific broker
        module for which the query is being updated.
        """
        self.object = self.get_object()
        return get_service_class(self.object.broker).form

    def get_form(self):
        """
        Returns an instance of the form to be used in this view.

        :returns: Form instance
        :rtype: django.forms.Form
        """
        form = super().get_form()
        form.helper.form_action = reverse(
            'public-update', kwargs={'pk': self.object.id}
        )
        return form

    def get_initial(self):
        """
        Returns the initial data to use for forms on this view. Initial data for this form consists of the name of
        the broker that the query is for.

        :returns: dict of initial values
        :rtype: dict
        """
        initial = super().get_initial()
        initial.update(self.object.parameters_as_dict)
        initial['broker'] = self.object.broker
        return initial

    def form_valid(self, form):
        """
        Saves the associated ``BrokerQuery`` and redirects to the ``BrokerQuery`` list.
        """
        form.save(query_id=self.object.id)
        return redirect(reverse('tom_alerts:list'))


class CustomPublicUserCreationForm(UserCreationForm):
    email = forms.EmailField(required=True)
    groups = forms.ModelMultipleChoiceField(Group.objects.all().exclude(name__in=['Public','AMON members']),
                                            required=True, widget=forms.CheckboxSelectMultiple,
                                            help_text='If your collaboration is an AMON member, use its account to access to private alerts.',
                                            )

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'groups')
        field_classes = {'username': UsernameField}

    def save(self, commit=True):
        user = super(forms.ModelForm, self).save(commit=False)
        if self.cleaned_data['password1']:
            user.set_password(self.cleaned_data["password1"])
        if commit:
            user.save()
            self.save_m2m()

        return user

    def __init__(self, *args, **kwargs):
        super(CustomPublicUserCreationForm, self).__init__(*args, **kwargs)
        self.fields["groups"].initial = (
            Group.objects.all().exclude(name__in=['Public','AMON members']).values_list(
                'id', flat=True
            )
        )

class PublicUserCreateView(CreateView):
    """
    View that handles ``User`` creation. Only for accessing public alerts.
    """
    template_name = 'tom_common/create_public_user.html'
    success_url = reverse_lazy('user-list')
    form_class = CustomPublicUserCreationForm

    def form_valid(self, form):
        """
        Called after form is validated. Creates the ``User`` and adds them to the public ``Group``.

        :param form: User creation form
        :type form: django.forms.Form
        """
        super().form_valid(form)
        group, _ = Group.objects.get_or_create(name='Public')
        group.user_set.add(self.object)
        group.save()
        return redirect(self.get_success_url())

class TargetCreateView(LoginRequiredMixin, CreateView):
    """
    View for creating a Target. Requires authentication.
    """

    model = Target
    fields = '__all__'

    def get_default_target_type(self):
        """
        Returns the user-configured target type specified in ``settings.py``, if it exists, otherwise returns sidereal

        :returns: User-configured target type or global default
        :rtype: str
        """
        try:
            return settings.TARGET_TYPE
        except AttributeError:
            return Target.SIDEREAL

    def get_target_type(self):
        """
        Gets the type of the target to be created from the query parameters. If none exists, use the default target
        type specified in ``settings.py``.

        :returns: target type
        :rtype: str
        """
        obj = self.request.GET or self.request.POST
        target_type = obj.get('type')
        # If None or some invalid value, use default target type
        if target_type not in (Target.SIDEREAL, Target.NON_SIDEREAL):
            target_type = self.get_default_target_type()
        return target_type

    def get_initial(self):
        """
        Returns the initial data to use for forms on this view.

        :returns: Dictionary with the following keys:

                  `type`: ``str``: Type of the target to be created

                  `groups`: ``QuerySet<Group>`` Groups available to the current user

        :rtype: dict
        """
        return {
            'type': self.get_target_type(),
            'groups': self.request.user.groups.all(),
            **dict(self.request.GET.items())
        }

    def get_context_data(self, **kwargs):
        """
        Inserts certain form data into the context dict.

        :returns: Dictionary with the following keys:

                  `type_choices`: ``tuple``: Tuple of 2-tuples of strings containing available target types in the TOM

                  `extra_form`: ``FormSet``: Django formset with fields for arbitrary key/value pairs
        :rtype: dict
        """
        context = super(TargetCreateView, self).get_context_data(**kwargs)
        context['type_choices'] = Target.TARGET_TYPES
        context['names_form'] = TargetNamesFormset(initial=[{'name': new_name}
                                                            for new_name
                                                            in self.request.GET.get('names', '').split(',')])
        context['extra_form'] = TargetExtraFormset()
        return context

    def get_form_class(self):
        """
        Return the form class to use in this view.

        :returns: form class for target creation
        :rtype: subclass of TargetCreateForm
        """
        target_type = self.get_target_type()
        self.initial['type'] = target_type
        if target_type == Target.SIDEREAL:
            return SiderealTargetCreateForm
        else:
            return NonSiderealTargetCreateForm

    def form_valid(self, form):
        """
        Runs after form validation. Creates the ``Target``, and creates any ``TargetName`` or ``TargetExtra`` objects,
        then runs the ``target_post_save`` hook and redirects to the success URL.

        :param form: Form data for target creation
        :type form: subclass of TargetCreateForm
        """
        super().form_valid(form)
        extra = TargetExtraFormset(self.request.POST)
        names = TargetNamesFormset(self.request.POST)
        if extra.is_valid() and names.is_valid():
            extra.instance = self.object
            extra.save()
            names.instance = self.object
            names.save()
        else:
            form.add_error(None, extra.errors)
            form.add_error(None, extra.non_form_errors())
            form.add_error(None, names.errors)
            form.add_error(None, names.non_form_errors())
            return super().form_invalid(form)
        logger.info('Target post save hook: %s created: %s', self.object, True)
        run_hook('target_post_save', target=self.object, created=True)
        return redirect(self.get_success_url())

    def get_form(self, *args, **kwargs):
        """
        Gets an instance of the ``TargetCreateForm`` and populates it with the groups available to the current user.

        :returns: instance of creation form
        :rtype: subclass of TargetCreateForm
        """
        form = super().get_form(*args, **kwargs)
        if self.request.user.is_superuser:
            form.fields['groups'].queryset = Group.objects.all()
        else:
            exclude_groups = ['Public']
            form.fields['groups'].queryset = self.request.user.groups.all().exclude(name__in=exclude_groups)
        return form


class TargetUpdateView(PermissionRequiredMixin, UpdateView):
    """
    View that handles updating a target. Requires authorization.
    """
    permission_required = 'tom_targets.change_target'
    model = Target
    fields = '__all__'

    def get_context_data(self, **kwargs):
        """
        Adds formset for ``TargetName`` and ``TargetExtra`` to the context.

        :returns: context object
        :rtype: dict
        """
        extra_field_names = [extra['name'] for extra in settings.EXTRA_FIELDS]
        context = super().get_context_data(**kwargs)
        context['names_form'] = TargetNamesFormset(instance=self.object)
        context['extra_form'] = TargetExtraFormset(
            instance=self.object,
            queryset=self.object.targetextra_set.exclude(key__in=extra_field_names)
        )
        return context

    @transaction.atomic
    def form_valid(self, form):
        """
        Runs after form validation. Validates and saves the ``TargetExtra`` and ``TargetName`` formsets, then calls the
        superclass implementation of ``form_valid``, which saves the ``Target``. If any forms are invalid, rolls back
        the changes.

        Saving is done in this order to ensure that new names/extras are available in the ``target_post_save`` hook.

        :param form: Form data for target update
        :type form: subclass of TargetCreateForm
        """
        extra = TargetExtraFormset(self.request.POST, instance=self.object)
        names = TargetNamesFormset(self.request.POST, instance=self.object)
        if extra.is_valid() and names.is_valid():
            extra.save()
            names.save()
        else:
            form.add_error(None, extra.errors)
            form.add_error(None, extra.non_form_errors())
            form.add_error(None, names.errors)
            form.add_error(None, names.non_form_errors())
            return super().form_invalid(form)
        super().form_valid(form)
        return redirect(self.get_success_url())

    def get_queryset(self, *args, **kwargs):
        """
        Returns the queryset that will be used to look up the Target by limiting the result to targets that the user is
        authorized to modify.

        :returns: Set of targets
        :rtype: QuerySet
        """
        return get_objects_for_user(self.request.user, 'tom_targets.change_target')

    def get_form_class(self):
        """
        Return the form class to use in this view.

        :returns: form class for target update
        :rtype: subclass of TargetCreateForm
        """
        if self.object.type == Target.SIDEREAL:
            return SiderealTargetCreateForm
        elif self.object.type == Target.NON_SIDEREAL:
            return NonSiderealTargetCreateForm

    def get_initial(self):
        """
        Returns the initial data to use for forms on this view. For the ``TargetUpdateView``, adds the groups that the
        target is a member of.

        :returns:
        :rtype: dict
        """
        initial = super().get_initial()
        initial['groups'] = get_groups_with_perms(self.get_object())
        return initial

    def get_form(self, *args, **kwargs):
        """
        Gets an instance of the ``TargetCreateForm`` and populates it with the groups available to the current user.

        :returns: instance of creation form
        :rtype: subclass of TargetCreateForm
        """
        form = super().get_form(*args, **kwargs)
        if self.request.user.is_superuser:
            form.fields['groups'].queryset = Group.objects.all()
        else:
            form.fields.pop('groups')
        return form

class CreateTargetFromAlertView(LoginRequiredMixin, View):
    """
    View that handles the creation of ``Target`` objects from a ``BrokerQuery`` result. Requires authentication.
    """

    def post(self, request, *args, **kwargs):
        """
        Handles the POST requests to this view. Creates a ``Target`` for each alert sent in the POST. Redirects to the
        ``TargetListView`` if multiple targets were created, and the ``TargetUpdateView`` if only one was created.
        Redirects to the ``RunQueryView`` if no ``Target`` objects. were successfully created.
        """
        query_id = self.request.POST['query_id']
        broker_name = self.request.POST['broker']
        broker_class = get_service_class(broker_name)
        alerts = self.request.POST.getlist('alerts')
        errors = []
        if not alerts:
            messages.warning(request, 'Please select at least one alert from which to create a target.')
            return redirect(reverse('tom_alerts:run', kwargs={'pk': query_id}))
        for alert_id in alerts:
            cached_alert = cache.get('alert_{}'.format(alert_id))
            if not cached_alert:
                messages.error(request, 'Could not create targets. Try re running the query again.')
                return redirect(reverse('tom_alerts:run', kwargs={'pk': query_id}))
            generic_alert = broker_class().to_generic_alert(json.loads(cached_alert))
            target = generic_alert.to_target()
            try:
                target.save()
                broker_class().process_reduced_data(target, json.loads(cached_alert))
                exclude_groups = ['Public']
                if broker_name in settings.PRIVATE_BROKERS:
                    exclude_groups += ['Public alerts']
                for group in request.user.groups.all().exclude(name__in=exclude_groups):
                    logger.info(group)
                    assign_perm('tom_targets.view_target', group, target)
                    assign_perm('tom_targets.change_target', group, target)
                    assign_perm('tom_targets.delete_target', group, target)
            except IntegrityError:
                messages.warning(request, f'Unable to save {target.name}, target with that name already exists.')
                errors.append(target.name)
        if (len(alerts) == len(errors)):
            return redirect(reverse('tom_alerts:run', kwargs={'pk': query_id}))
        elif (len(alerts) == 1):
            return redirect(reverse(
                'tom_targets:update', kwargs={'pk': target.id})
            )
        else:
            return redirect(reverse(
                'tom_targets:list')
            )

class AmonTargetDetailView(Raise403PermissionRequiredMixin, DetailView):
    """
    View that handles the display of the target details. Requires authorization.
    """
    permission_required = 'tom_targets.view_target'
    model = Target

    def get_context_data(self, *args, **kwargs):
        """
        Adds the ``DataProductUploadForm`` to the context and prepopulates the hidden fields.

        :returns: context object
        :rtype: dict
        """
        context = super().get_context_data(*args, **kwargs)
        #data_product_upload_form = DataProductUploadForm(
        #    initial={
        #        'target': self.get_object(),
        #        'referrer': reverse('targets:detail', args=(self.get_object().id,))
        #    }
        #)
        #context['data_product_form'] = data_product_upload_form
        observing_strategy_form = RunStrategyForm(initial={'target': self.get_object()})
        if any(self.request.GET.get(x) for x in ['observing_strategy', 'cadence_strategy', 'cadence_frequency']):
            initial = {'target': self.object}
            initial.update(self.request.GET)
            observing_strategy_form = RunStrategyForm(
                initial=initial
            )
        observing_strategy_form.fields['target'].widget = HiddenInput()
        context['observing_strategy_form'] = observing_strategy_form
        return context

    def get(self, request, *args, **kwargs):
        """
        Handles the GET requests to this view. If update_status is passed into the query parameters, calls the
        updatestatus management command to query for new statuses for ``ObservationRecord`` objects associated with this
        target.

        :param request: the request object passed to this view
        :type request: HTTPRequest
        """
        update_status = request.GET.get('update_status', False)
        if update_status:
            if not request.user.is_authenticated:
                return redirect(reverse('login'))
            target_id = kwargs.get('pk', None)
            out = StringIO()
            call_command('updatestatus', target_id=target_id, stdout=out)
            messages.info(request, out.getvalue())
            add_hint(request, mark_safe(
                              'Did you know updating observation statuses can be automated? Learn how in'
                              '<a href=https://tom-toolkit.readthedocs.io/en/stable/customization/automation.html>'
                              ' the docs.</a>'))
            #return redirect(reverse('targets:detail', args=(target_id,)))
            return redirect(reverse('tom_targets:detail', args=(target_id,)))

        run_strategy_form = RunStrategyForm(request.GET)
        if run_strategy_form.is_valid():
            obs_strat = ObservingStrategy.objects.get(pk=run_strategy_form.cleaned_data['observing_strategy'].id)
            target_id = kwargs.get('pk', None)
            params = urlencode(obs_strat.parameters_as_dict)
            params += urlencode(request.GET)
            return redirect(
                reverse('tom_observations:create',
                        args=(obs_strat.facility,)) + f'?target_id={self.get_object().id}&' + params)

        return super().get(request, *args, **kwargs)
