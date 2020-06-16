"""django URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.urls import path, include, re_path
from django.views.generic import TemplateView
from django.views.defaults import page_not_found
import django
from .views import AboutView, PrivateAmonRunQueryView, PublicAmonRunQueryView, AmonBrokerQueryListView
from .views import PublicBrokerQueryCreateView, PublicBrokerQueryUpdateView, PublicBrokerQueryDeleteView, PublicUserCreateView
from .views import TargetCreateView, TargetUpdateView, CreateTargetFromAlertView, AmonTargetDetailView

urlpatterns = [
    # path('about/', TemplateView.as_view(template_name='about.html'), name='about'),
    path('about/', AboutView.as_view(), name='about'),
    path('alerts/query/<int:pk>/run/', PublicAmonRunQueryView.as_view(), name='PublicAmonRunQuery'),
    path('alerts/query/<int:pk>/run/private/', PrivateAmonRunQueryView.as_view(), name='PrivateAmonRunQuery'),
    path('alerts/query/list/', AmonBrokerQueryListView.as_view(), name='AmonBrokerQueryList'),
    re_path('alerts/query/.*/run/', django.views.defaults.page_not_found, {'exception': Exception('Not Found')}),
    path('alert/query/create/public/', PublicBrokerQueryCreateView.as_view(), name='public-create'),
    path('alert/query/<int:pk>/update/public/', PublicBrokerQueryUpdateView.as_view(), name='public-update'),
    path('alert/query/<int:pk>/delete/', PublicBrokerQueryDeleteView.as_view(), name='public-delete'),
    path('sign_up/', PublicUserCreateView.as_view(), name='public-user-create'),
    path('targets/create/', TargetCreateView.as_view(), name='create-target'),
    path('alerts/alert/create/', CreateTargetFromAlertView.as_view(), name='create-target-from-alert'),
    path('targets/<int:pk>/update/', TargetUpdateView.as_view(), name='update-target'),
    path('targets/<int:pk>/', AmonTargetDetailView.as_view(), name='amon-target-detail'),
    path('', include('tom_common.urls')), # NB This should be last
]
# NB to change a page, add url here before `path('', include('tom_common.urls'))`, change the name to the one you wrote here in the html that call this url (like a buton)
