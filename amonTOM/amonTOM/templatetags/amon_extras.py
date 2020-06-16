import json

from django import template
from django.conf import settings
from django.contrib.auth.models import Group
from django.shortcuts import reverse
from guardian.shortcuts import get_groups_with_perms

from tom_dataproducts.forms import DataProductUploadForm
from tom_observations.models import ObservationRecord
from tom_targets.models import Target

from astropy import units as u
from django.conf import settings
from plotly import offline
from plotly import graph_objs as go


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
            form.fields['groups'].queryset = get_groups_with_perms(obj)
    return {'data_product_form': form}

@register.inclusion_tag('tom_targets/partials/target_distribution.html')
def target_distribution_amon(targets):
    """
    Displays a plot showing on a map the locations of all sidereal targets in the TOM.
    """
    locations = targets.filter(type=Target.SIDEREAL).values_list('ra', 'dec', 'name')
    data = [
        dict(
            lon=[l[0] for l in locations],
            lat=[l[1] for l in locations],
            text=[l[2] for l in locations],
            hoverinfo='lon+lat+text',
            mode='markers',
            type='scattergeo'
        ),
        dict(
            lon=list(range(0, 360, 60))+[180]*4,
            lat=[0]*6+[-60, -30, 30, 60],
            text=list(range(0, 360, 60))+[-60, -30, 30, 60],
            hoverinfo='none',
            mode='text',
            type='scattergeo'
        )
    ]
    layout = {
        'title': 'Target Distribution (sidereal)',
        'hovermode': 'closest',
        'showlegend': False,
        'height': 350,
        'margin': {"r":0,"t":25,"l":0,"b":0},
        'geo': {
            'projection': {
                'type': 'mollweide',
            },
            'showcoastlines': False,
            'showland': False,
            'lonaxis': {
                'showgrid': True,
                'range': [0, 360],
            },
            'lataxis': {
                'showgrid': True,
                'range': [-90, 90],
            },
        }
    }
    figure = offline.plot(go.Figure(data=data, layout=layout), output_type='div', show_link=False)
    return {'figure': figure}

@register.inclusion_tag('tom_observations/partials/observation_distribution.html')
def observation_distribution_amon(observations):
    """
    Displays a plot showing on a map the locations of all observations recorded in the TOM.
    """

    # "distinct" query is not supported, must manually find distinct observation per target
    sorted_observations = observations.order_by('scheduled_end')  # ascending so that only the max is preserved
    observation_targets = {}
    for obs in sorted_observations:
        observation_targets[obs.target_id] = (obs.status, obs.terminal)

    observation_no_status = [t for t in observation_targets.keys()
                             if not observation_targets[t][0]]  # status==""
    observation_terminal = [t for t in observation_targets.keys()
                            if observation_targets[t][0]
                            and observation_targets[t][1]]  # status!="" and terminal
    observation_non_terminal = [t for t in observation_targets.keys()
                                if observation_targets[t][0]
                                and not observation_targets[t][1]]  # status!="" and not terminal

    targets_no_status = Target.objects.filter(pk__in=observation_no_status)
    targets_terminal = Target.objects.filter(pk__in=observation_terminal)
    targets_non_terminal = Target.objects.filter(pk__in=observation_non_terminal)

    locations_no_status = targets_no_status.filter(type=Target.SIDEREAL).values_list('ra', 'dec', 'name')
    locations_terminal = targets_terminal.filter(type=Target.SIDEREAL).values_list('ra', 'dec', 'name')
    locations_non_terminal = targets_non_terminal.filter(type=Target.SIDEREAL).values_list('ra', 'dec', 'name')

    data = [
        dict(
            lon=[l[0] for l in locations_no_status],
            lat=[l[1] for l in locations_no_status],
            text=[l[2] for l in locations_no_status],
            hoverinfo='lon+lat+text',
            mode='markers',
            marker=dict(color='rgba(90, 90, 90, .8)'),
            type='scattergeo'
        ),
        dict(
            lon=[l[0] for l in locations_non_terminal],
            lat=[l[1] for l in locations_non_terminal],
            text=[l[2] for l in locations_non_terminal],
            hoverinfo='lon+lat+text',
            mode='markers',
            marker=dict(color='rgba(152, 0, 0, .8)'),
            type='scattergeo'
        ),
        dict(
            lon=[l[0] for l in locations_terminal],
            lat=[l[1] for l in locations_terminal],
            text=[l[2] for l in locations_terminal],
            hoverinfo='lon+lat+text',
            mode='markers',
            marker=dict(color='rgba(0, 152, 0, .8)'),
            type='scattergeo'
        ),
        dict(
            lon=list(range(0, 360, 60))+[180]*4,
            lat=[0]*6+[-60, -30, 30, 60],
            text=list(range(0, 360, 60))+[-60, -30, 30, 60],
            hoverinfo='none',
            mode='text',
            type='scattergeo'
        )
    ]
    layout = {
        'title': 'Observation Distribution (sidereal)',
        'hovermode': 'closest',
        'showlegend': False,
        'height': 350,
        'margin': {"r":0,"t":25,"l":0,"b":0},
        'geo': {
            'projection': {
                'type': 'mollweide',
            },
            'showcoastlines': False,
            'showland': False,
            'lonaxis': {
                'showgrid': True,
                'range': [0, 360],
            },
            'lataxis': {
                'showgrid': True,
                'range': [-90, 90],
            },
        }
    }
    figure = offline.plot(go.Figure(data=data, layout=layout), output_type='div', show_link=False)
    return {'figure': figure}
