import logging

from tom_alerts.alerts import GenericQueryForm, GenericAlert, GenericBroker
from tom_alerts.models import BrokerQuery
from tom_targets.models import Target
from crispy_forms.layout import Layout, Div, Fieldset, HTML
from dataclasses import dataclass
import datetime
import pytz
from astropy.time import Time

from django import forms
from django.conf import settings
import requests

from astropy import units as u
from astropy.coordinates import SkyCoord

from sshtunnel import SSHTunnelForwarder
import MySQLdb as db
import pandas as pd
import numpy as np
from numpy import sin, cos
import os

logger = logging.getLogger(__name__)

# Example of broker: https://github.com/TOMToolkit/tom_base/blob/master/tom_alerts/brokers/mars.py

# ssh variables
#host = '3.13.26.235' # This is AWS dev
host = '3.132.123.216'
localhost = '127.0.0.1'
ssh_username = 'ubuntu'
ssh_private_key = settings.SSH_KEY_FILE

# database variables
user='amon'
password=settings.AMON_DB_PASSWORD
#database='AMON_test' # for AWS dev
database='AMON_test2'

def equatorial_to_galactic(ra, dec):
    direction = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    return direction.galactic.l.degree, direction.galactic.b.degree

def J2000(time,ra,dec,degrees=True):
    """
        Convert celestial coordinates in current epoch to J2000 epoch.
        Receives time in UTC, and cooridnates in degrees as default
    """
    try:
        t = Time(time,format='isot',scale='utc')
    except ValueError:
        t = Time(time,format='iso',scale='utc')
    JD = t.jd #Julian day
    JC = (JD - 2451545.0)/36525.0 #Julian Century

    zeta   = (0.6406161 + (0.0000839 + 0.0000050*JC)*JC)*JC
    z      = (0.6406161 + (0.0003041 + 0.0000051*JC)*JC)*JC
    theta  = (0.5567530 - (0.0001185 + 0.0000116*JC)*JC)*JC

    zeta = np.deg2rad(zeta)
    z = np.deg2rad(z)
    theta = np.deg2rad(theta)

    cosTheta= np.cos(theta)
    sinTheta = np.sin(theta)

    if degrees:
        ra = np.deg2rad(ra)
        dec = np.deg2rad(dec)

    sinDec = np.sin(dec)
    cosDec = np.cos(dec)
    cosRA = np.cos(ra-z)

    dec2000 = -cosRA * sinTheta * cosDec + cosTheta * sinDec
    dec2000 = np.arcsin(dec2000)
    dec2000 = np.rad2deg(dec2000)

    ra2000 = 0.0
    if dec2000 == 90.0:
        return ra2000, dec2000

    s1 = np.sin(ra-z)*cosDec
    c2 = cosRA * cosTheta * cosDec + sinTheta * sinDec
    ra2000 = np.arctan2(s1,c2) - zeta
    ra2000 = np.rad2deg(ra2000)

    if ra2000<0:
        ra2000 += 360.

    return ra2000, dec2000

def query(q):
    """Connect to the AMON machine and query the mysql database"""
    with SSHTunnelForwarder(
          (host, 22),
          ssh_username=ssh_username,
          ssh_private_key=ssh_private_key,
          # ssh_password=ssh_password,
          remote_bind_address=(localhost, 3306)
     ) as server:
        for attempt in range(3):
            try:
                conn = db.connect(host=localhost,
                    port=server.local_bind_port,
                    user=user,
                    passwd=password,
                    db=database,
                    connect_timeout=30)
            except:
                if attempt <= 1:
                    logger.info('Connection error, trying again.')
                else:
                    logger.info('Connection error. Stopping.')
                continue
            else:
                logger.info('Connection succeded.')
                break
        return pd.read_sql_query(q, conn)

@dataclass
class ICGoldBronzeAlert:
    """
    dataclass representing an alert in order to display it in the UI.
    """

    timestamp: datetime.datetime
    id: int
    name: str
    ra: float
    dec: float
    l: float
    b: float
    # mag: float
    # score: float
    url: str
    energy: float
    charge: float
    signalness: float
    far: float
    src_error: float
    src_error90: float
    stream: int
    rev: int

    def to_target(self):
        """
        Returns a Target instance for an object defined by an alert.

        :returns: representation of object for an alert
        :rtype: `Target`
        """
        return Target(
            name=self.name,
            type='SIDEREAL',
            ra=self.ra,
            dec=self.dec
        )

class ICGoldBronzeBrokerForm(GenericQueryForm):
    """Set up the fields of the broker form"""
    evt_num = forms.IntegerField(required=False)
    stream_list = ((24, 'Gold'), (25, 'Bronze'), (10, 'HESE'), (11, 'EHE'))
    # stream = forms.ChoiceField(choices=streams,
    #     required=False, label='Streams'
    # )
    streams = forms.MultipleChoiceField(choices=stream_list, required=True,
        widget=forms.CheckboxSelectMultiple, label='Streams',
        help_text="HESE and EHE ended in June 2019 when Gold and Bronze started.",
    )
    time__gt = forms.CharField(required=False, label='Time Lower Bound',
        widget=forms.TextInput(attrs={'type': 'date'}),
    )
    time__lt = forms.CharField(required=False, label='Time Upper Bound',
        widget=forms.TextInput(attrs={'type': 'date'})
    )
    time__since = forms.IntegerField(required=False, label='Time Since',
        help_text='Alerts younger than this number of seconds',
    )
    jd__gt = forms.FloatField(required=False, label='JD Lower Bound')
    jd__lt = forms.FloatField(required=False, label='JD Upper Bound')
    cone = forms.CharField(required=False, label='Cone Search',
        help_text='RA,Dec,radius in degrees'
    )
    # objectcone = forms.CharField(required=False, label='Object Cone Search',
    #     help_text='Object name,radius in degrees'
    # )
    ra__gt = forms.FloatField(required=False, label='RA [°] Lower Bound', help_text='[0,360]')
    ra__lt = forms.FloatField(required=False, label='RA [°] Upper Bound', help_text='[0,360]')
    dec__gt = forms.FloatField(required=False, label='Dec [°] Lower Bound', help_text='[-90,90]')
    dec__lt = forms.FloatField(required=False, label='Dec [°] Upper Bound', help_text='[-90,90]')
    l__gt = forms.FloatField(required=False, label='l [°] Lower Bound', help_text='[0,360]')
    l__lt = forms.FloatField(required=False, label='l [°] Upper Bound', help_text='[0,360]')
    b__gt = forms.FloatField(required=False, label='b [°] Lower Bound', help_text='[-90,90]')
    b__lt = forms.FloatField(required=False, label='b [°] Upper Bound', help_text='[-90,90]')
    err__lt = forms.FloatField(required=False, label='90% Uncertainty Upper Bound')
    err__gt = forms.FloatField(required=False, label='90% Uncertainty Lower Bound')
    err50__lt = forms.FloatField(required=False, label='50% Uncertainty Upper Bound')
    err50__gt = forms.FloatField(required=False, label='50% Uncertainty Lower Bound')
    ener__lt = forms.FloatField(required=False, label='Energy [GeV] Upper Bound')
    ener__gt = forms.FloatField(required=False, label='Energy [GeV] Lower Bound')
    qtot__lt = forms.FloatField(required=False, label='Charge [pe] Upper Bound', help_text='HESE or EHE')
    qtot__gt = forms.FloatField(required=False, label='Charge [pe] Lower Bound', help_text='HESE or EHE')
    sig__lt = forms.FloatField(required=False, label='Signalness [dn] Upper Bound')
    sig__gt = forms.FloatField(required=False, label='Signalness [dn] Lower Bound')
    far__lt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Upper Bound')
    far__gt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Lower Bound')


    
    def __init__(self, *args, **kwargs):
        """Structure the form page"""
        super().__init__(*args, **kwargs)
        self.helper.layout = Layout(
            HTML('''
                <p>
                Please see the <a href="https://gcn.gsfc.nasa.gov/amon.html">GCN AMON documentation</a>
                for a detailed description of the <a href="https://gcn.gsfc.nasa.gov/doc/IceCube_High_Energy_Neutrino_Track_Alerts_v2.pdf">Gold and Bronze</a>, <a href="https://gcn.gsfc.nasa.gov/doc/Public_Doc_AMON_IceCube_GCN_Alerts_Oct2016_v7.pdf">HESE</a> or <a href="https://gcn.gsfc.nasa.gov/doc/AMON_IceCube_EHE_alerts_Oct31_2016.pdf">EHE</a> alerts.
                </p>
            '''),
            self.common_layout,
            'evt_num', 'streams',
            Fieldset('Time based filters', 'time__since',
                Div(
                    Div('time__gt', 'jd__gt', css_class='col',),
                    Div('time__lt', 'jd__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset(
                'Energy based filters',
                Div(
                    Div('ener__gt', css_class='col',),
                    Div('ener__lt', css_class='col',),
                    css_class="form-row",
                ),
                Div(
                    Div('qtot__gt', css_class='col',),
                    Div('qtot__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset(
                'Signalness based filters',
                Div(
                    Div('sig__gt', 'far__gt', css_class='col',),
                    Div('sig__lt', 'far__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset('Location based filters', 'cone',
                # 'objectcone',
                Div(
                    Div('ra__gt', 'dec__gt', 'l__gt', 'b__gt', css_class='col',),
                    Div('ra__lt', 'dec__lt', 'l__lt', 'b__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset(
                'Angular uncertainty based filters',
                Div(
                    Div('err__gt', 'err50__gt', css_class='col',),
                    Div('err__lt', 'err50__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
        )


class ICGoldBronzeBroker(GenericBroker):
    """Make the mysql selection string from the form output and query the database"""
    name = 'IceCube Track'
    form = ICGoldBronzeBrokerForm

    @classmethod
    def fetch_alerts(clazz, params_sel): # NB should not give access to subthreshold events. check that alert and event tables are ok
        
        # List of param from "event" table of DB
        param_evttab = [['evt_num', 'id='],
                       ['time__since', 'time>'],
                       ['time__gt', 'time>'],
                       ['time__lt', 'time<'],
                       ['jd__gt', 'time>'],
                       ['jd__lt', 'time<'],
                       ['cone', ''],
                       # ['objectcone', ''],
                       ['ra__gt', 'RA>'],
                       ['ra__lt', 'RA<'],
                       ['dec__gt', '`Dec`>'],
                       ['dec__lt', '`Dec`<'],
                       ['l__gt', '>'],
                       ['l__lt', '<'],
                       ['b__gt', '>'],
                       ['b__lt', '<'],
                      ]

        # List of param from "parameter" table of DB
        param_partab = [['ener__gt', 'energy', 'gt'],
                       ['ener__lt', 'energy', 'lt'],
                       ['qtot__gt', 'causalqtot', 'gt'],
                       ['qtot__lt', 'causalqtot', 'lt'],
                       ['qtot__gt', 'qtot', 'gt'],
                       ['qtot__lt', 'qtot', 'lt'],
                       ['sig__gt', 'signalness', 'gt'],
                       ['sig__lt', 'signalness', 'lt'],
                       ['sig__gt', 'signal_trackness', 'gt'],
                       ['sig__lt', 'signal_trackness', 'lt'],
                       ['far__gt', 'far', 'gt'],
                       ['far__lt', 'far', 'lt'],
                       ['err__gt', 'src_error90', 'gt'],
                       ['err__lt', 'src_error90', 'lt'],
                       ['err50__gt', 'src_error', 'gt'],
                       ['err50__lt', 'src_error', 'lt'],
                      ]

        # Streams
        streams = params_sel['streams']
        selection = "SELECT * FROM event WHERE type='observation' AND"
        selection += " ("
        for count, stream in enumerate(streams):
            if stream == '24' or stream == '25':
                selection += "(time>'2019-06-19' AND eventStreamConfig_stream={stream})".format(stream=stream)
            elif stream == '10':
                selection += "(time>'2016-04-26' AND eventStreamConfig_stream={stream})".format(stream=stream)
            elif stream == '11':
                selection += "(time>'2016-07-30' AND eventStreamConfig_stream={stream})".format(stream=stream)
            else:
                selection +="eventStreamConfig_stream={stream}".format(stream=stream)
            if count < len(streams)-1:
                selection += " OR "
        selection += ")"
        # print(selection, file=open('/var/www/amonTOM/amonTOM/print.txt', 'a'))

        # Other params_sel
        for param_tom, db_condit in param_evttab:
            if params_sel[param_tom] is None or params_sel[param_tom] is '':
                continue
            if param_tom == 'time__since':
                now = datetime.datetime.now(pytz.UTC)
                timedelta = datetime.timedelta(seconds=params_sel[param_tom])
                time_since = now - timedelta
                selection += " AND {db_condit}'{time_since}'".format(db_condit=db_condit, time_since=str(time_since))
            elif param_tom in ['jd__gt', 'jd__lt']:
                time = Time(params_sel[param_tom], format='mjd')
                timedate = time.isot
                selection += " AND {db_condit}'{timedate}'".format(db_condit=db_condit, timedate=timedate)
            elif param_tom in ['b__gt', 'b__lt']:
                # NB longitude and latitude in db are not set (always 0 and -90)
                # ...so I need to write the conversion when querying the database (and I do the longitude selection later)
                Dec_NGP = str(np.deg2rad(27.12825))
                Ra_NGP = str(np.deg2rad(192.85948))
                # Lon_NCP = str(np.deg2rad(122.93192))
                conv = 'SIN({dec_NGP})*SIN(RADIANS(`Dec`))+COS({dec_NGP})*COS(RADIANS(`Dec`))*COS(RADIANS(`RA`)-{ra_NGP})'.format(dec_NGP=Dec_NGP, ra_NGP=Ra_NGP)
                selection += " AND {conv}{db_condit}SIN(RADIANS('{param_tom}'))".format(conv=conv, db_condit=db_condit, param_tom=params_sel[param_tom])
            elif param_tom == 'cone':
                center_ra = float(params_sel[param_tom].split(',')[0])
                center_dec = float(params_sel[param_tom].split(',')[1])
                cone_radius = float(params_sel[param_tom].split(',')[2])
                selection += " AND SQRT(POW(`RA`-{ra_center},2) + POW(`Dec`-{dec_center},2))<{radius}".format(ra_center=center_ra, dec_center=center_dec, radius=cone_radius)
            elif params_sel[param_tom] is not '' and param_tom not in ['l__gt', 'l__lt']: # in general
                selection += " AND {db_condit}'{param_tom}'".format(db_condit=db_condit, param_tom=params_sel[param_tom]) # '' are necessary for the time
        if '10' in streams and '11' in streams:
            selection += " AND NOT (id=1282906888376 AND eventStreamConfig_stream=11)" # NB Hardcoded: There is one event which is HESE and EHE at the same time, so I select it as HESE only and modify it after
        selection += ";"
        
        logger.info(selection)
        df = query(selection)
        df['time'] = df['time'].astype(str)
        df['energy'] = None # empty 
        df['signalness'] = None # empty 
        df['signal_trackness'] = None # empty 
        df['far'] = None # empty 
        df['src_error'] = None # empty 
        df['src_error90'] = None # empty 
        df['causalqtot'] = None # empty 
        df['qtot'] = None # empty 
        df['stream'] = ""
        df['stream'][df['eventStreamConfig_stream'] == 24] = "Gold"
        df['stream'][df['eventStreamConfig_stream'] == 25] = "Bronze"
        df['stream'][df['eventStreamConfig_stream'] == 10] = "HESE"
        df['stream'][df['eventStreamConfig_stream'] == 11] = "EHE"
        df['stream'][df['id']==1282906888376] = "HESE EHE"

        dic_alerts = df.to_dict('records')

        for param_tom in ['l__gt', 'l__lt']:
            # NB longitude and latitude in db are not set (always 0 and -90)
            # ...so I do the selection afterwards for longitude
            if params_sel[param_tom] is None:
                continue
            ra_list = [dic_alerts[i]['RA'] for i in range(len(dic_alerts))]
            dec_list = [dic_alerts[i]['Dec'] for i in range(len(dic_alerts))]
            l, b = equatorial_to_galactic(ra_list, dec_list)
            if '__gt' in param_tom:
                not_sel = l > params_sel[param_tom]
            else:
                not_sel = l < params_sel[param_tom]
            dic_alerts = np.array(dic_alerts)[not_sel]

        ids = [alert['id'] for alert in dic_alerts]
        ids_repeated = list(set([id for id in ids if ids.count(id) > 1]))
        index_del = []
        for id in ids_repeated:
            indices_thisid = np.where(np.array(ids) == id)[0]
            revs = []
            for index in indices_thisid:
                revs += [dic_alerts[index]['rev']]
            for index, rev in zip(indices_thisid, revs):
                if rev < max(revs):
                    index_del += [index]
        dic_alerts = np.delete(dic_alerts, index_del)


        # Query the "parameter" table to get energy, far, error etc
        n_match = len(dic_alerts)
        if n_match > 0:
            selection_param = "SELECT * FROM parameter WHERE "
            selection_param += " ("
            # for count, stream in enumerate(streams):
            #     selection_param +=" event_eventStreamConfig_stream={stream}".format(stream=stream)
            #     if count < len(streams)-1:
            #         selection_param += " OR "
            # selection_param += ") AND ("
            for index, alert in enumerate(dic_alerts):
                if index!=0:
                    selection_param += " OR "
                # selection_param += "event_id={id}".format(id=alert['id'])
                selection_param += "(event_id={id} AND event_eventStreamConfig_stream={stream})".format(id=alert['id'], stream=alert['eventStreamConfig_stream'])

            begin = True
            for param_tom, param_name, db_condit in param_partab:
                if begin:
                    selection_param += ") AND ("
                else:
                    selection_param += " OR "
                selection_param += "name='{param_name}'".format(param_name=param_name)
                begin = False
            selection_param += ");"

            logger.info(selection_param)
            df = query(selection_param)
            # df['time'] = df['time'].astype(str)
            dic_params = df.to_dict('records')
            # logger.info(dic_params)
            # Get event id of params not following conditions to remove them from dic_alerts
            rej_ids = []
            for param in dic_params:
                if param['event_id'] in rej_ids:
                    continue
                for param_tom, param_name, db_condit in param_partab:
                    if (params_sel[param_tom] is not None and params_sel[param_tom] is not '' # if there is a condition on this param
                            and param['name'] == param_name and
                            ((param['value'] < params_sel[param_tom] and db_condit == 'gt') # if not follow gt condition
                            or (param['value'] > params_sel[param_tom] and db_condit == 'lt'))): # or if not follow lt condition
                        rej_ids += [param['event_id']]
            index_del = []
            for index, alert in enumerate(dic_alerts):
                if alert['id'] in rej_ids:
                    index_del += [index]
            dic_alerts = np.delete(dic_alerts, index_del)

            # Put param values (energy, far, etc) in dic_alerts
            ids = np.array([alert['id'] for alert in dic_alerts])
            revs = np.array([alert['rev'] for alert in dic_alerts])
            for param in dic_params:
                rev_sel = (ids == param['event_id']) * (revs == param['event_rev'])
                if not np.any(rev_sel):
                    continue
                index = np.where(rev_sel)[0][0]
                if param['name'] != 'event_id':
                    if param['name'] == 'energy':
                        dic_alerts[index][param['name']]=param['value']/1000. # To get TeV
                    else:
                        dic_alerts[index][param['name']]=param['value']

        return iter([alert for alert in dic_alerts])

    @classmethod
    def to_generic_alert(clazz, alert):
        isHESE = ("HESE" in alert['stream'])
        isGoldBronze = (alert['stream'] == "Gold") or (alert['stream'] == "Bronze")
        url = "https://gcn.gsfc.nasa.gov/notices_amon{_g_b}/{id1}_{id2}.amon".format(
                        _g_b='_g_b' if isGoldBronze else '',
                        id1=str(alert['id'])[0:6] if isGoldBronze else str(alert['id'])[6:], # run_num if GB evt_num else
                        id2=str(alert['id'])[6:] if isGoldBronze else str(alert['id'])[0:6], # evt_num if GB run_num else
                        )
        return ICGoldBronzeAlert(
            timestamp=alert['time'],
            url=url,
            id=alert['id'],
            name="IC_"+str(alert['id']),
            ra=alert['RA'],
            dec=alert['Dec'],
            l=equatorial_to_galactic(alert['RA'], alert['Dec'])[0],
            b=equatorial_to_galactic(alert['RA'], alert['Dec'])[1],
            # mag=None,
            # score=None,
            energy=alert['energy'],
            charge=alert['causalqtot'] if isHESE else alert['qtot'],
            signalness=alert['signal_trackness'] if isHESE else alert['signalness'],
            far=alert['far'],
            src_error=alert['src_error'],
            src_error90=alert['src_error90'] if isHESE or isGoldBronze else None,
            stream=alert['stream'],
            rev=alert['rev'],
        )


@dataclass
class ICCascadeAlert:
    """
    dataclass representing an alert in order to display it in the UI.
    """

    timestamp: datetime.datetime
    id: int
    name: str
    ra: float
    dec: float
    l: float
    b: float
    # mag: float
    # score: float
    url: str
    energy: float
    #charge: float
    signalness: float
    far: float
    src_error: float
    src_error90: float
    stream: int
    fits_url: str
    png_url: str
    rev: int

    def to_target(self):
        """
        Returns a Target instance for an object defined by an alert.

        :returns: representation of object for an alert
        :rtype: `Target`
        """
        return Target(
            name=self.name,
            type='SIDEREAL',
            ra=self.ra,
            dec=self.dec
        )

class ICCascadeBrokerForm(GenericQueryForm):
    """Set up the fields of the broker form"""
    evt_num = forms.IntegerField(required=False)
    time__gt = forms.CharField(required=False, label='Time Lower Bound',
        widget=forms.TextInput(attrs={'type': 'date'}),
    )
    time__lt = forms.CharField(required=False, label='Time Upper Bound',
        widget=forms.TextInput(attrs={'type': 'date'})
    )
    time__since = forms.IntegerField(required=False, label='Time Since',
        help_text='Alerts younger than this number of seconds',
    )
    jd__gt = forms.FloatField(required=False, label='JD Lower Bound')
    jd__lt = forms.FloatField(required=False, label='JD Upper Bound')
    cone = forms.CharField(required=False, label='Cone Search',
        help_text='RA,Dec,radius in degrees'
    )
    ra__gt = forms.FloatField(required=False, label='RA [°] Lower Bound', help_text='[0,360]')
    ra__lt = forms.FloatField(required=False, label='RA [°] Upper Bound', help_text='[0,360]')
    dec__gt = forms.FloatField(required=False, label='Dec [°] Lower Bound', help_text='[-90,90]')
    dec__lt = forms.FloatField(required=False, label='Dec [°] Upper Bound', help_text='[-90,90]')
    l__gt = forms.FloatField(required=False, label='l [°] Lower Bound', help_text='[0,360]')
    l__lt = forms.FloatField(required=False, label='l [°] Upper Bound', help_text='[0,360]')
    b__gt = forms.FloatField(required=False, label='b [°] Lower Bound', help_text='[-90,90]')
    b__lt = forms.FloatField(required=False, label='b [°] Upper Bound', help_text='[-90,90]')
    err__lt = forms.FloatField(required=False, label='90% Uncertainty Upper Bound')
    err__gt = forms.FloatField(required=False, label='90% Uncertainty Lower Bound')
    err50__lt = forms.FloatField(required=False, label='50% Uncertainty Upper Bound')
    err50__gt = forms.FloatField(required=False, label='50% Uncertainty Lower Bound')
    ener__lt = forms.FloatField(required=False, label='Energy [GeV] Upper Bound')
    ener__gt = forms.FloatField(required=False, label='Energy [GeV] Lower Bound')
    sig__lt = forms.FloatField(required=False, label='Signalness [dn] Upper Bound')
    sig__gt = forms.FloatField(required=False, label='Signalness [dn] Lower Bound')
    far__lt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Upper Bound')
    far__gt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Lower Bound')


    
    def __init__(self, *args, **kwargs):
        """Structure the form page"""
        super().__init__(*args, **kwargs)
        self.helper.layout = Layout(
            HTML('''
                <p>
                Please see the <a href="https://gcn.gsfc.nasa.gov/doc/High_Energy_Neutrino_Cascade_Alerts.pdf">Cascade alerts documentation</a>
                for a detailed description of <a href="https://gcn.gsfc.nasa.gov/amon.html">the stream</a>.
                </p>
            '''),
            self.common_layout,
            'evt_num', 'streams',
            Fieldset('Time based filters', 'time__since',
                Div(
                    Div('time__gt', 'jd__gt', css_class='col',),
                    Div('time__lt', 'jd__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset(
                'Energy based filters',
                Div(
                    Div('ener__gt', css_class='col',),
                    Div('ener__lt', css_class='col',),
                    css_class="form-row",
                ),
            ),
            Fieldset(
                'Signalness based filters',
                Div(
                    Div('sig__gt', 'far__gt', css_class='col',),
                    Div('sig__lt', 'far__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset('Location based filters', 'cone',
                # 'objectcone',
                Div(
                    Div('ra__gt', 'dec__gt', 'l__gt', 'b__gt', css_class='col',),
                    Div('ra__lt', 'dec__lt', 'l__lt', 'b__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset(
                'Angular uncertainty based filters',
                Div(
                    Div('err__gt', 'err50__gt', css_class='col',),
                    Div('err__lt', 'err50__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
        )


class ICCascadeBroker(GenericBroker):
    """Make the mysql selection string from the form output and query the database"""
    name = 'IceCube Cascade'
    form = ICCascadeBrokerForm

    @classmethod
    def fetch_alerts(clazz, params_sel): # NB should not give access to subthreshold events. check that alert and event tables are ok
        
        # List of param from "event" table of DB
        param_evttab = [['evt_num', 'id='],
                       ['time__since', 'time>'],
                       ['time__gt', 'time>'],
                       ['time__lt', 'time<'],
                       ['jd__gt', 'time>'],
                       ['jd__lt', 'time<'],
                       ['cone', ''],
                       ['ra__gt', 'RA>'],
                       ['ra__lt', 'RA<'],
                       ['dec__gt', '`Dec`>'],
                       ['dec__lt', '`Dec`<'],
                       ['l__gt', '>'],
                       ['l__lt', '<'],
                       ['b__gt', '>'],
                       ['b__lt', '<'],
                      ]

        # List of param from "parameter" table of DB
        param_partab = [['ener__gt', 'energy', 'gt'],
                       ['ener__lt', 'energy', 'lt'],
                       ['sig__gt', 'signalness', 'gt'],
                       ['sig__lt', 'signalness', 'lt'],
                       ['far__gt', 'far', 'gt'],
                       ['far__lt', 'far', 'lt'],
                       ['err__gt', 'src_error90', 'gt'],
                       ['err__lt', 'src_error90', 'lt'],
                       ['err50__gt', 'src_error', 'gt'],
                       ['err50__lt', 'src_error', 'lt'],
                      ]

        # Streams
        #streams = params_sel['streams']
        stream = 26 #str(streams[0])
        selection = "SELECT * FROM event WHERE type='observation' AND"
        selection += " eventStreamConfig_stream={stream}".format(stream=stream)
        #selection += " AND signalness!=-1"

        # Other params_sel
        for param_tom, db_condit in param_evttab:
            if params_sel[param_tom] is None or params_sel[param_tom] is '':
                continue
            if param_tom == 'time__since':
                now = datetime.datetime.now(pytz.UTC)
                timedelta = datetime.timedelta(seconds=params_sel[param_tom])
                time_since = now - timedelta
                selection += " AND {db_condit}'{time_since}'".format(db_condit=db_condit, time_since=str(time_since))
            elif param_tom in ['jd__gt', 'jd__lt']:
                time = Time(params_sel[param_tom], format='mjd')
                timedate = time.isot
                selection += " AND {db_condit}'{timedate}'".format(db_condit=db_condit, timedate=timedate)
            elif param_tom in ['b__gt', 'b__lt']:
                # NB longitude and latitude in db are not set (always 0 and -90)
                # ...so I need to write the conversion when querying the database (and I do the longitude selection later)
                Dec_NGP = str(np.deg2rad(27.12825))
                Ra_NGP = str(np.deg2rad(192.85948))
                # Lon_NCP = str(np.deg2rad(122.93192))
                conv = 'SIN({dec_NGP})*SIN(RADIANS(`Dec`))+COS({dec_NGP})*COS(RADIANS(`Dec`))*COS(RADIANS(`RA`)-{ra_NGP})'.format(dec_NGP=Dec_NGP, ra_NGP=Ra_NGP)
                selection += " AND {conv}{db_condit}SIN(RADIANS('{param_tom}'))".format(conv=conv, db_condit=db_condit, param_tom=params_sel[param_tom])
            elif param_tom == 'cone':
                center_ra = float(params_sel[param_tom].split(',')[0])
                center_dec = float(params_sel[param_tom].split(',')[1])
                cone_radius = float(params_sel[param_tom].split(',')[2])
                selection += " AND SQRT(POW(`RA`-{ra_center},2) + POW(`Dec`-{dec_center},2))<{radius}".format(ra_center=center_ra, dec_center=center_dec, radius=cone_radius)
            elif params_sel[param_tom] is not '' and param_tom not in ['l__gt', 'l__lt']: # in general
                selection += " AND {db_condit}'{param_tom}'".format(db_condit=db_condit, param_tom=params_sel[param_tom]) # '' are necessary for the time
        selection += ";"
        
        logger.info(selection)
        df = query(selection)
        df['time'] = df['time'].astype(str)
        df['energy'] = None # empty 
        df['signalness'] = None # empty 
        df['far'] = None # empty 
        df['src_error'] = None # empty 
        df['src_error90'] = None # empty 
        df['name'] = None # empty
        df['fits_url'] = None # empty
        df['png_url'] = None # empty
        df['stream'] = "Cascade"

        dic_alerts = df.to_dict('records')

        for param_tom in ['l__gt', 'l__lt']:
            # NB longitude and latitude in db are not set (always 0 and -90)
            # ...so I do the selection afterwards for longitude
            if params_sel[param_tom] is None:
                continue
            ra_list = [dic_alerts[i]['RA'] for i in range(len(dic_alerts))]
            dec_list = [dic_alerts[i]['Dec'] for i in range(len(dic_alerts))]
            l, b = equatorial_to_galactic(ra_list, dec_list)
            if '__gt' in param_tom:
                not_sel = l > params_sel[param_tom]
            else:
                not_sel = l < params_sel[param_tom]
            dic_alerts = np.array(dic_alerts)[not_sel]

        ids = [alert['id'] for alert in dic_alerts]
        ids_repeated = list(set([id for id in ids if ids.count(id) > 1]))
        index_del = []
        for id in ids_repeated:
            indices_thisid = np.where(np.array(ids) == id)[0]
            revs = []
            for index in indices_thisid:
                revs += [dic_alerts[index]['rev']]
            for index, rev in zip(indices_thisid, revs):
                if rev < max(revs):
                    index_del += [index]
        dic_alerts = np.delete(dic_alerts, index_del)


        n_match = len(dic_alerts)
        if n_match > 0:
            # Query the "parameter" table to get energy, far, error etc
            selection_param = "SELECT * FROM parameter WHERE "
            selection_param += " ("
            for index, alert in enumerate(dic_alerts):
                if index!=0:
                    selection_param += " OR "
                selection_param += "(event_id={id} AND event_eventStreamConfig_stream={stream})".format(id=alert['id'], stream=alert['eventStreamConfig_stream'])

            begin = True
            for param_tom, param_name, db_condit in param_partab:
                if begin:
                    selection_param += ") AND ("
                else:
                    selection_param += " OR "
                selection_param += "name='{param_name}'".format(param_name=param_name)
                begin = False
            selection_param += " OR name like 'IceCubeCascade-%'"

            selection_param += ");"

            logger.info(selection_param)
            df = query(selection_param)
            dic_params = df.to_dict('records')
            # Get event id of params not following conditions to remove them from dic_alerts
            rej_ids = []
            for param in dic_params:
                if param['event_id'] in rej_ids:
                    continue
                if param['name'] == 'signalness' and param['value'] < 0: # rej events with sig = -1 ie subthreshold events
                    rej_ids += [param['event_id']]
                for param_tom, param_name, db_condit in param_partab:
                    if (params_sel[param_tom] is not None and params_sel[param_tom] is not '' # if there is a condition on this param
                            and param['name'] == param_name and
                            ((param['value'] < params_sel[param_tom] and db_condit == 'gt') # if not follow gt condition
                            or (param['value'] > params_sel[param_tom] and db_condit == 'lt'))): # or if not follow lt condition
                        rej_ids += [param['event_id']]
                        
            index_del = []
            for index, alert in enumerate(dic_alerts):
                if alert['id'] in rej_ids:
                    index_del += [index]
            dic_alerts = np.delete(dic_alerts, index_del)

            # Put param values (energy, far, etc) in dic_alerts
            ids = np.array([alert['id'] for alert in dic_alerts])
            revs = np.array([alert['rev'] for alert in dic_alerts])
            logger.info(dic_params)
            for param in dic_params:
                rev_sel = (ids == param['event_id']) * (revs == param['event_rev'])
                if not np.any(rev_sel):
                    continue
                index = np.where(rev_sel)[0][0]
                if param['name'] != 'event_id':
                    if param['name'] == 'energy':
                        dic_alerts[index][param['name']]=param['value']/1000. # To get TeV
                    elif param['name'][:15] == 'IceCubeCascade-':
                        dic_alerts[index]['name']=param['name']
                    else:
                        dic_alerts[index][param['name']]=param['value']

            # Query the "skyMapEvent" table to get fits and png skymap url
            selection_skymap = "SELECT * FROM skyMapEvent WHERE"
            selection_skymap += " ("
            for index, alert in enumerate(dic_alerts):
                if index!=0:
                    selection_skymap += " OR "
                selection_skymap += "(event_id={id} AND event_eventStreamConfig_stream={stream})".format(id=alert['id'], stream=alert['eventStreamConfig_stream'])

            selection_skymap += ");"

            if len(dic_alerts) == 0: # when there is no events selected
                selection_skymap = "SELECT * FROM skyMapEvent limit 0"

            logger.info(selection_skymap)
            df = query(selection_skymap)
            dic_skymaps = df.to_dict('records')
            # Get event id of params not following conditions to remove them from dic_alerts
            rej_ids = []
            for skymap in dic_skymaps:
                if skymap['event_id'] in rej_ids:
                    continue

            # Put skymap urls in dic_alerts
            ids = np.array([alert['id'] for alert in dic_alerts])
            revs = np.array([alert['rev'] for alert in dic_alerts])
            for skymap in dic_skymaps:
                rev_sel = (ids == skymap['event_id']) * (revs == skymap['event_rev'])
                if not np.any(rev_sel):
                    continue
                index = np.where(rev_sel)[0][0]
                if skymap['location'][-4:] == '.png':
                    dic_alerts[index]['png_url']=skymap['location']
                elif skymap['location'][-5:] == '.fits':
                    dic_alerts[index]['fits_url']=skymap['location']

        return iter([alert for alert in dic_alerts])

    @classmethod
    def to_generic_alert(clazz, alert):
        return ICCascadeAlert(
            timestamp=alert['time'],
            url="https://gcn.gsfc.nasa.gov/notices_amon_icecube_cascade/{run_num}_{evt_num}.amon".format(run_num=str(alert['id'])[0:6], evt_num=str(int(str(alert['id'])[6:]))),

            id=alert['id'],
            name=alert['name'],
            ra=alert['RA'],
            dec=alert['Dec'],
            l=equatorial_to_galactic(alert['RA'], alert['Dec'])[0],
            b=equatorial_to_galactic(alert['RA'], alert['Dec'])[1],
            energy=alert['energy'],
            signalness=alert['signalness'],
            far=alert['far'],
            src_error=alert['src_error'],
            src_error90=alert['src_error90'],
            stream=alert['stream'],
            rev=alert['rev'],
            fits_url=alert['fits_url'],
            png_url=alert['png_url'],
        )




@dataclass
class NuEMAlert:
    """
    dataclass representing an alert in order to display it in the UI.
    """

    timestamp: datetime.datetime
    id: int
    name: str
    ra: float
    dec: float
    l: float
    b: float
    src_error90: float
    nevents: int
    url: str
    far: float
    delta_t: float
    sigma_t: float
    stream: str
    rev: int

    def to_target(self):
        """
        Returns a Target instance for an object defined by an alert.

        :returns: representation of object for an alert
        :rtype: `Target`
        """
        return Target(
            name=self.name,
            type='SIDEREAL',
            ra=self.ra,
            dec=self.dec
        )


class NuEMAlertBrokerForm(GenericQueryForm):
    """Set up the fields of the broker form"""
    evt_num = forms.IntegerField(required=False)
    stream_list = ((1, 'IceCube-HAWC'), (8, 'ANTARES-Fermi'))
    # stream = forms.ChoiceField(choices=streams,
    #     required=False, label='Streams'
    # )
    streams = forms.MultipleChoiceField(choices=stream_list, required=True, widget=forms.CheckboxSelectMultiple, label='Streams')
    time__gt = forms.CharField(required=False, label='Time Lower Bound',
        widget=forms.TextInput(attrs={'type': 'date'})
    )
    time__lt = forms.CharField(required=False, label='Time Upper Bound',
        widget=forms.TextInput(attrs={'type': 'date'})
    )
    time__since = forms.IntegerField(required=False, label='Time Since',
        help_text='Alerts younger than this number of seconds'
    )
    jd__gt = forms.FloatField(required=False, label='JD Lower Bound')
    jd__lt = forms.FloatField(required=False, label='JD Upper Bound')
    cone = forms.CharField(required=False, label='Cone Search',
        help_text='RA,Dec,radius in degrees'
    )
    # objectcone = forms.CharField(required=False, label='Object Cone Search',
    #     help_text='Object name,radius in degrees'
    # )
    ra__gt = forms.FloatField(required=False, label='RA [°] Lower Bound', help_text='[0,360]')
    ra__lt = forms.FloatField(required=False, label='RA [°] Upper Bound', help_text='[0,360]')
    dec__gt = forms.FloatField(required=False, label='Dec [°] Lower Bound', help_text='[-90,90]')
    dec__lt = forms.FloatField(required=False, label='Dec [°] Upper Bound', help_text='[-90,90]')
    l__gt = forms.FloatField(required=False, label='l [°] Lower Bound', help_text='[0,360]')
    l__lt = forms.FloatField(required=False, label='l [°] Upper Bound', help_text='[0,360]')
    b__gt = forms.FloatField(required=False, label='b [°] Lower Bound', help_text='[-90,90]')
    b__lt = forms.FloatField(required=False, label='b [°] Upper Bound', help_text='[-90,90]')
    sig__lt = forms.FloatField(required=False, label='Signalness [dn] Upper Bound')
    sig__gt = forms.FloatField(required=False, label='Signalness [dn] Lower Bound')
    far__lt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Upper Bound')
    far__gt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Lower Bound')


    
    def __init__(self, *args, **kwargs):
        """Structure the form page"""
        super().__init__(*args, **kwargs)
        self.helper.layout = Layout(
            HTML('''
                <p>
                Please see the <a href="https://gcn.gsfc.nasa.gov/doc/gamma_nu.pdf">gamma-neutrino coincidence alerts documentation</a>
                for a detailed description of <a href="https://gcn.gsfc.nasa.gov/amon.html">the stream</a>.
                </p>
            '''),
            self.common_layout,
            'evt_num', 'streams',
            Fieldset('Time based filters', 'time__since',
                Div(
                    Div('time__gt', 'jd__gt', css_class='col',),
                    Div('time__lt', 'jd__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            # Fieldset(
            #     'Energy based filters',
            #     Div(
            #         Div('ener__gt', css_class='col',),
            #         Div('ener__lt', css_class='col',),
            #         css_class="form-row",
            #     )
            # ),
            Fieldset(
                'Signalness based filters',
                Div(
                    Div('sig__gt', 'far__gt', css_class='col',),
                    Div('sig__lt', 'far__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset('Location based filters', 'cone',
                # 'objectcone',
                Div(
                    Div('ra__gt', 'dec__gt', 'l__gt', 'b__gt', css_class='col',),
                    Div('ra__lt', 'dec__lt', 'l__lt', 'b__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            # Fieldset(
            #     'Angular uncertainty based filters',
            #     Div(
            #         Div('err__gt', 'err50__gt', css_class='col',),
            #         Div('err__lt', 'err50__lt', css_class='col',),
            #         css_class="form-row",
            #     )
            # ),
        )


class NuEMAlertBroker(GenericBroker):
    """Make the mysql selection string from the form output and query the database"""
    name = 'Nu-EM'
    form = NuEMAlertBrokerForm

    @classmethod
    def fetch_alerts(clazz, params_sel): # NB should not give access to subthreshold events. check that alert and event tables are ok
        
        # List of param from "event" table of DB
        param_alrtab = [['evt_num', 'id='],
                       ['time__since', 'time>'],
                       ['time__gt', 'time>'],
                       ['time__lt', 'time<'],
                       ['jd__gt', 'time>'],
                       ['jd__lt', 'time<'],
                       ['sig__gt', 'pvalue>'],
                       ['sig__lt', 'pvalue<'],
                       ['far__gt', 'false_pos>'],
                       ['far__lt', 'false_pos<'],
                       ['cone', ''],
                       # ['objectcone', ''],
                       ['ra__gt', 'RA>'],
                       ['ra__lt', 'RA<'],
                       ['dec__gt', '`Dec`>'],
                       ['dec__lt', '`Dec`<'],
                       ['l__gt', '>'],
                       ['l__lt', '<'],
                       ['b__gt', '>'],
                       ['b__lt', '<'],
                      ]

        
        # alert table contains only streams IC-HAWC=1 and Antares-Fermi=8
        selection = "SELECT * FROM alert WHERE type='observation' AND"
        
        streams = params_sel['streams']
        stream_cond = {'1': " AND false_pos<=4", '8': " AND ((time>'2019-05-01' AND time<'2020-07-07' AND false_pos<=4/(3600*24*365.25)) OR (false_pos<=4 AND time>'2020-07-07'))"} # stream: threshold. Before july 2020 ant-fermi false_pos was in s-1, now in yr-1
        selection += " ("
        for count, stream in enumerate(streams):
            selection += "(alertConfig_stream={stream} {stream_cond})".format(stream=stream, stream_cond=stream_cond[stream])
            if count < len(streams)-1:
                selection += " OR "
        selection += ")"
        # print(selection, file=open('/var/www/amonTOM/amonTOM/print.txt', 'a'))

        # Other params_sel
        for param_tom, db_condit in param_alrtab:
            if params_sel[param_tom] is None or params_sel[param_tom] is '':
                continue
            if param_tom == 'time__since':
                now = datetime.datetime.now(pytz.UTC)
                timedelta = datetime.timedelta(seconds=params_sel[param_tom])
                time_since = now - timedelta
                selection += " AND {db_condit}'{time_since}'".format(db_condit=db_condit, time_since=str(time_since))
            elif param_tom in ['jd__gt', 'jd__lt']:
                time = Time(params_sel[param_tom], format='mjd')
                timedate = time.isot
                selection += " AND {db_condit}'{timedate}'".format(db_condit=db_condit, timedate=timedate)
            elif param_tom in ['b__gt', 'b__lt']:
                # NB longitude and latitude in db are not set (always 0 and -90)
                # ...so I need to write the conversion when querying the database (and I do the longitude selection later)
                Dec_NGP = str(np.deg2rad(27.12825))
                Ra_NGP = str(np.deg2rad(192.85948))
                # Lon_NCP = str(np.deg2rad(122.93192))
                conv = 'SIN({dec_NGP})*SIN(RADIANS(`Dec`))+COS({dec_NGP})*COS(RADIANS(`Dec`))*COS(RADIANS(`RA`)-{ra_NGP})'.format(dec_NGP=Dec_NGP, ra_NGP=Ra_NGP)
                selection += " AND {conv}{db_condit}SIN(RADIANS('{param_tom}'))".format(conv=conv, db_condit=db_condit, param_tom=params_sel[param_tom])
            elif param_tom == 'cone':
                center_ra = float(params_sel[param_tom].split(',')[0])
                center_dec = float(params_sel[param_tom].split(',')[1])
                cone_radius = float(params_sel[param_tom].split(',')[2])
                selection += " AND SQRT(POW(`RA`-{ra_center},2) + POW(`Dec`-{dec_center},2))<{radius}".format(ra_center=center_ra, dec_center=center_dec, radius=cone_radius)
            elif params_sel[param_tom] is not '' and param_tom not in ['l__gt', 'l__lt']:
                selection += " AND {db_condit}'{param_tom}'".format(db_condit=db_condit, param_tom=params_sel[param_tom]) # '' are necessary for the time
        
        selection += ";"

        df = query(selection)
        #df['time'] = df['time'].astype(str)
        df['sigmaT'][df['alertConfig_stream'] == 1] = ' ' # empty
        select_false_pos_sec = (df['alertConfig_stream'] == 8) * (df['time']<datetime.datetime(2020, 7, 7))
        df['false_pos'][select_false_pos_sec] = df['false_pos'][select_false_pos_sec]*3600.*24.*365.25
        df['sigmaR'][df['alertConfig_stream'] == 8] = df['sigmaR'][df['alertConfig_stream'] == 8]*2.146 # NB Colin may want to use something more conservative (see slack conv we had on 25 Mar 2020)

        alert_config = query("SELECT * FROM alertConfig WHERE (stream=1 OR stream=8)") # Used to know when run number change (which correspond to a change in the analysis)
        df['run'] = np.zeros(len(df['sigmaT']))
        for index in range(len(alert_config['stream'])):
            df['run'][np.array(df['alertConfig_stream']==alert_config['stream'][index]) * np.array(df['time'].to_numpy()>alert_config['validStart'][index].to_numpy())] = alert_config['rev'][index]

        df['alertConfig_stream'][df['alertConfig_stream'] == 1] = "IC-HAWC"
        df['alertConfig_stream'][df['alertConfig_stream'] == 8] = "ANT-Fermi"
        df['time'] = df['time'].astype(str)

        dic_alerts = df.to_dict('records')

        for param_tom in ['l__gt', 'l__lt']:
            # NB longitude and latitude in db are not set (always 0 and -90)
            # ...so I do the selection afterwards for longitude
            if params_sel[param_tom] is None:
                continue
            ra_list = [dic_alerts[i]['RA'] for i in range(len(dic_alerts))]
            dec_list = [dic_alerts[i]['Dec'] for i in range(len(dic_alerts))]
            l, b = equatorial_to_galactic(ra_list, dec_list)
            if '__gt' in param_tom:
                not_sel = l > params_sel[param_tom]
            else:
                not_sel = l < params_sel[param_tom]
            dic_alerts = np.array(dic_alerts)[not_sel]
        
        # Use last rev
        ids = [alert['id'] for alert in dic_alerts]
        ids_repeated = list(set([id for id in ids if ids.count(id) > 1]))
        index_del = []
        for id in ids_repeated:
            indices_thisid = np.where(np.array(ids) == id)[0]
            revs = []
            for index in indices_thisid:
                revs += [dic_alerts[index]['rev']]
            for index, rev in zip(indices_thisid, revs):
                if rev < max(revs):
                    index_del += [index]
        dic_alerts = np.delete(dic_alerts, index_del)

        return iter([alert for alert in dic_alerts])

    @classmethod
    def to_generic_alert(clazz, alert):
        url = "https://gcn.gsfc.nasa.gov/notices_amon_nu_em/{evt}_{run}.amon".format(
                evt=str(alert['id']), run=str(int(alert['run']))) # NB run number is 0 until we change something to the analysis. If that is the case I should change it here.
        circ_first_alerts = { # there was circulars before notices were ready
                37610: 'https://gcn.gsfc.nasa.gov/gcn3/26963.gcn3',
                1436255238: 'https://gcn.gsfc.nasa.gov/gcn3/26005.gcn3',
                1440223121: 'https://gcn.gsfc.nasa.gov/gcn3/26674.gcn3',
                1441204786: 'https://gcn.gsfc.nasa.gov/gcn3/26915.gcn3',
                }
        if alert['id'] in circ_first_alerts.keys():
            url = circ_first_alerts[alert['id']]
        return NuEMAlert(
            timestamp=alert['time'],
            url=url,
            id=alert['id'],
            name="NuEM_"+str(alert['id']),
            ra=alert['RA'],
            dec=alert['Dec'],
            l=equatorial_to_galactic(alert['RA'], alert['Dec'])[0],
            b=equatorial_to_galactic(alert['RA'], alert['Dec'])[1],
            src_error90=alert['sigmaR'],
            nevents=alert['nevents'],
            far=alert['false_pos'],
            delta_t=alert['deltaT'],
            sigma_t=alert['sigmaT'],
            stream=alert['alertConfig_stream'],
            rev=alert['rev'],
        )


@dataclass
class HAWCGRBAlert:
    """
    dataclass representing an alert in order to display it in the UI.
    """

    timestamp: datetime.datetime
    id: int
    name: str
    ra: float
    dec: float
    l: float
    b: float
    # mag: float
    # score: float
    url: str
    deltaT: float
    pvalue: float
    far: float
    src_error90: float
    stream: int
    rev: int

    def to_target(self):
        """
        Returns a Target instance for an object defined by an alert.

        :returns: representation of object for an alert
        :rtype: `Target`
        """
        return Target(
            name=self.name,
            type='SIDEREAL',
            ra=self.ra,
            dec=self.dec
        )

class HAWCGRBBrokerForm(GenericQueryForm):
    """Set up the fields of the broker form"""
    evt_num = forms.IntegerField(required=False)
    # streams = ((8, 'HAWC GRB'))
    # stream = forms.ChoiceField(choices=streams,
    #     required=False, label='Streams'
    # )
    time__gt = forms.CharField(required=False, label='Time Lower Bound',
        widget=forms.TextInput(attrs={'type': 'date'})
    )
    time__lt = forms.CharField(required=False, label='Time Upper Bound',
        widget=forms.TextInput(attrs={'type': 'date'})
    )
    time__since = forms.IntegerField(required=False, label='Time Since',
        help_text='Alerts younger than this number of seconds',
    )
    jd__gt = forms.FloatField(required=False, label='JD Lower Bound')
    jd__lt = forms.FloatField(required=False, label='JD Upper Bound')
    cone = forms.CharField(required=False, label='Cone Search',
        help_text='RA,Dec,radius in degrees'
    )
    # objectcone = forms.CharField(required=False, label='Object Cone Search',
    #     help_text='Object name,radius in degrees'
    # )
    ra__gt = forms.FloatField(required=False, label='RA [°] Lower Bound', help_text='[0,360]')
    ra__lt = forms.FloatField(required=False, label='RA [°] Upper Bound', help_text='[0,360]')
    dec__gt = forms.FloatField(required=False, label='Dec [°] Lower Bound', help_text='[-90,90]')
    dec__lt = forms.FloatField(required=False, label='Dec [°] Upper Bound', help_text='[-90,90]')
    l__gt = forms.FloatField(required=False, label='l [°] Lower Bound', help_text='[0,360]')
    l__lt = forms.FloatField(required=False, label='l [°] Upper Bound', help_text='[0,360]')
    b__gt = forms.FloatField(required=False, label='b [°] Lower Bound', help_text='[-90,90]')
    b__lt = forms.FloatField(required=False, label='b [°] Upper Bound', help_text='[-90,90]')
    err__lt = forms.FloatField(required=False, label='90% Uncertainty Upper Bound')
    err__gt = forms.FloatField(required=False, label='90% Uncertainty Lower Bound')
    # err50__lt = forms.FloatField(required=False, label='50% Uncertainty Upper Bound')
    # err50__gt = forms.FloatField(required=False, label='50% Uncertainty Lower Bound')
    # ener__lt = forms.FloatField(required=False, label='Energy [GeV] Upper Bound')
    # ener__gt = forms.FloatField(required=False, label='Energy [GeV] Lower Bound')
    pval__lt = forms.FloatField(required=False, label='Pvalue Upper Bound', help_text='[0,1]')
    pval__gt = forms.FloatField(required=False, label='Pvalue Lower Bound', help_text='[0,1]')
    far__lt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Upper Bound')
    far__gt = forms.FloatField(required=False, label=r'False Alarm Rate [yr⁻¹] Lower Bound')


    def __init__(self, *args, **kwargs):
        """Structure the form page"""
        super().__init__(*args, **kwargs)
        self.helper.layout = Layout(
            HTML('''
                <p>
                Please see the <a href="https://gcn.gsfc.nasa.gov/amon.html">GCN AMON documentation</a>,
                in particular <a href="https://gcn.gsfc.nasa.gov/doc/hawc_grb_alerts.pdf">here</a>, for a detailed description of these alerts.
                </p>
            '''),
            self.common_layout,
            'evt_num', 'stream',
            Fieldset('Time based filters', 'time__since',
                Div(
                    Div('time__gt', 'jd__gt', css_class='col',),
                    Div('time__lt', 'jd__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            # Fieldset(
            #     'Energy based filters',
            #     Div(
            #         Div('ener__gt', css_class='col',),
            #         Div('ener__lt', css_class='col',),
            #         css_class="form-row",
            #     )
            # ),
            Fieldset(
                'Signalness based filters',
                Div(
                    Div('pval__gt', 'far__gt', css_class='col',),
                    Div('pval__lt', 'far__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset('Location based filters', 'cone',
                # 'objectcone',
                Div(
                    Div('ra__gt', 'dec__gt', 'l__gt', 'b__gt', css_class='col',),
                    Div('ra__lt', 'dec__lt', 'l__lt', 'b__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
            Fieldset(
                'Angular uncertainty based filters',
                Div(
                    Div('err__gt', css_class='col',),
                    Div('err__lt', css_class='col',),
                    css_class="form-row",
                )
            ),
        )


class HAWCGRBBroker(GenericBroker):
    """Make the mysql selection string from the form output and query the database"""
    name = 'HAWC-GRB'
    form = HAWCGRBBrokerForm

    @classmethod
    def fetch_alerts(clazz, params_sel): # NB should not give access to subthreshold events. check that alert and event tables are ok
        
        # List of param from "event" table of DB
        param_evttab = [['evt_num', 'id='],
                       ['time__since', 'time>'],
                       ['time__gt', 'time>'],
                       ['time__lt', 'time<'],
                       ['jd__gt', 'time>'],
                       ['jd__lt', 'time<'],
                       ['cone', ''],
                       # ['objectcone', ''],
                       ['ra__gt', 'RA>'],
                       ['ra__lt', 'RA<'],
                       ['dec__gt', '`Dec`>'],
                       ['dec__lt', '`Dec`<'],
                       ['l__gt', '>'],
                       ['l__lt', '<'],
                       ['b__gt', '>'],
                       ['b__lt', '<'],
                       ['pval__gt', 'pvalue>'],
                       ['pval__lt', 'pvalue<'],
                       ['far__gt', 'false_pos>'],
                       ['far__lt', 'false_pos<'],
                      ]

        # Streams
        stream = 8 # HAWC GRB stream number
        selection = "SELECT * FROM event WHERE type='observation' AND"
        selection += " eventStreamConfig_stream={stream}".format(stream=stream)
        selection += " AND false_pos<=12" # NB this is to get only public events
        selection += " AND time>'2019-08-01'" # NB this is to get only public events

        # Other params_sel
        for param_tom, db_condit in param_evttab:
            if params_sel[param_tom] is None or params_sel[param_tom] is '':
                continue
            if param_tom == 'time__since':
                now = datetime.datetime.now(pytz.UTC)
                timedelta = datetime.timedelta(seconds=params_sel[param_tom])
                time_since = now - timedelta
                selection += " AND {db_condit}'{time_since}'".format(db_condit=db_condit, time_since=str(time_since))
            elif param_tom in ['jd__gt', 'jd__lt']:
                time = Time(params_sel[param_tom], format='mjd')
                timedate = time.isot
                selection += " AND {db_condit}'{timedate}'".format(db_condit=db_condit, timedate=timedate)
            elif param_tom in ['b__gt', 'b__lt']:
                # NB longitude and latitude in db are not set (always 0 and -90)
                # ...so I need to write the conversion when querying the database (and I do the longitude selection later)
                Dec_NGP = str(np.deg2rad(27.12825))
                Ra_NGP = str(np.deg2rad(192.85948))
                # Lon_NCP = str(np.deg2rad(122.93192))
                conv = 'SIN({dec_NGP})*SIN(RADIANS(`Dec`))+COS({dec_NGP})*COS(RADIANS(`Dec`))*COS(RADIANS(`RA`)-{ra_NGP})'.format(dec_NGP=Dec_NGP, ra_NGP=Ra_NGP)
                selection += " AND {conv}{db_condit}SIN(RADIANS('{param_tom}'))".format(conv=conv, db_condit=db_condit, param_tom=params_sel[param_tom])
            elif param_tom == 'cone':
                center_ra = float(params_sel[param_tom].split(',')[0])
                center_dec = float(params_sel[param_tom].split(',')[1])
                cone_radius = float(params_sel[param_tom].split(',')[2])
                selection += " AND SQRT(POW(`RA`-{ra_center},2) + POW(`Dec`-{dec_center},2))<{radius}".format(ra_center=center_ra, dec_center=center_dec, radius=cone_radius)
            elif params_sel[param_tom] is not '' and param_tom not in ['l__gt', 'l__lt']: # in general
                selection += " AND {db_condit}'{param_tom}'".format(db_condit=db_condit, param_tom=params_sel[param_tom]) # '' are necessary for the time
        selection += ";"
        
        logger.info(selection)
        df = query(selection)
        df['time'] = df['time'].astype(str)
        dic_alerts = df.to_dict('records')
        for i in range(len(dic_alerts)):
            dic_alerts[i]['RA'], dic_alerts[i]['Dec'] = J2000(dic_alerts[i]['time'], dic_alerts[i]['RA'], dic_alerts[i]['Dec'])

        for param_tom in ['l__gt', 'l__lt']:
            # NB longitude and latitude in db are not set (always 0 and -90)
            # ...so I do the selection afterwards for longitude
            if params_sel[param_tom] is None:
                continue
            ra_list = [dic_alerts[i]['RA'] for i in range(len(dic_alerts))]
            dec_list = [dic_alerts[i]['Dec'] for i in range(len(dic_alerts))]
            l, b = equatorial_to_galactic(ra_list, dec_list)
            if '__gt' in param_tom:
                not_sel = l > params_sel[param_tom]
            else:
                not_sel = l < params_sel[param_tom]
            dic_alerts = np.array(dic_alerts)[not_sel]

        # Get last rev only
        ids = [alert['id'] for alert in dic_alerts]
        ids_repeated = list(set([id for id in ids if ids.count(id) > 1]))
        index_del = []
        for id in ids_repeated:
            indices_thisid = np.where(np.array(ids) == id)[0]
            revs = []
            for index in indices_thisid:
                revs += [dic_alerts[index]['rev']]
            for index, rev in zip(indices_thisid, revs):
                if rev < max(revs):
                    index_del += [index]
        dic_alerts = np.delete(dic_alerts, index_del)

        return iter([alert for alert in dic_alerts])

    @classmethod
    def to_generic_alert(clazz, alert):
        return HAWCGRBAlert(
            timestamp=alert['time'],
            url="https://gcn.gsfc.nasa.gov/notices_amon_hawc/{run_num}_{evt_num}.amon".format(run_num=str(alert['id'])[0:-8], evt_num=str(int(str(alert['id'])[-8:]))),
            id=alert['id'],
            name="HawcGRB_"+str(alert['id']),
            ra=alert['RA'],
            dec=alert['Dec'],
            l=equatorial_to_galactic(alert['RA'], alert['Dec'])[0],
            b=equatorial_to_galactic(alert['RA'], alert['Dec'])[1],
            # mag=None,
            # score=None,
            deltaT=alert['deltaT'],
            pvalue=alert['pvalue'],
            far=alert['false_pos'],
            src_error90=alert['sigmaR'],
            stream=alert['eventStreamConfig_stream'],
            rev=alert['rev'],
        )
