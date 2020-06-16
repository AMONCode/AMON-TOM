import json

from django import template
from django.conf import settings
from django.contrib.auth.models import Group
from django.core.paginator import Paginator
from django.shortcuts import reverse
from datetime import datetime
from guardian.shortcuts import get_objects_for_user

from plotly import offline
import plotly.graph_objs as go

from tom_dataproducts.forms import DataProductUploadForm
from tom_dataproducts.models import DataProduct, ReducedDatum
from tom_dataproducts.processors.data_serializers import SpectrumSerializer
from tom_observations.models import ObservationRecord
from tom_targets.models import Target

register = template.Library()

@register.inclusion_tag('tom_dataproducts/partials/upload_dataproduct.html', takes_context=True)
def upload_dataproduct_amon(context, obj):
    user = context['user']
    initial = {}
    
    if isinstance(obj, Target):
        initial['target'] = obj
        initial['referrer'] = reverse('tom_targets:detail', args=(obj.id,))
    elif isinstance(obj, ObservationRecord):
        initial['observation_record'] = obj
        initial['referrer'] = reverse('tom_observations:detail', args=(obj.id,))
    form = DataProductUploadForm(initial=initial)
    if not settings.TARGET_PERMISSIONS_ONLY:
        if user.is_superuser:
            form.fields['groups'].queryset = Group.objects.all()
        else:
            form.fields['groups'].queryset = user.groups.all().exclude(name__in=['Public alerts'])
    return {'data_product_form': form}
