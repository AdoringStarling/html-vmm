import dash
import dash_bootstrap_components as dbc
from dash import html
from dash import dcc
import pandas as pd
import plotly.graph_objects as go
import numpy as np
import os
from geoseismo import *
from base64 import b64encode
import io

#Area del estudio
los=-72.72
loi=-75.41
lai=4.15
las=10.39

#Se cargan los datos de elevacion estos fueron descargados en https://portal.opentopography.org/datasets
# Global Bathymetry and Topography at 15 Arc Sec: SRTM15+ V2.1  
df_topo   =pd.read_csv('datasets\dem_srtm15arcs_area_LBG.xyz',delimiter=',',header=None,decimal='.')
df_topo   =df_topo[(df_topo[1]>lai)&(df_topo[1]<las)&(df_topo[0]>loi)&(df_topo[0]<los)] #Filtros previos
mesh_topo = (df_topo.pivot(index=1, columns=0,values=2))
z_topo,x_topo,y_topo=mesh_topo.values,mesh_topo.columns,mesh_topo.index

#Base de datos de sismos convertidos a csv desde http://bdrsnc.sgc.gov.co/paginas1/catalogo/Consulta_Valle_Medio/valle_medio.php
#df_sismos=pd.read_csv("datasets/reporte_1.csv")#,delimiter=';',decimal=',')
df_sismos=pd.read_csv(r'datasets\reporte_LBG.csv')
df_sismos['FECHA - HORA UTC']=df_sismos['Fecha  (UTC)'].astype(str)+' '+df_sismos['Hora  (UTC)'].astype(str)
df_sismos.rename(columns = {'Latitud(°)':'LATITUD (°)', 
                                'Longitud(°)':'LONGITUD (°)',
                                'Profundidad(Km)':'PROF. (Km)',
                                'Magnitud':'MAGNITUD',
                                'Tipo Magnitud':'TIPO MAGNITUD',
                                'Rms(Seg)':'RMS (Seg)',
                                'Gap(°)':'GAP (°)',
                                'Error  Latitud(Km)':'ERROR LATITUD (Km)',
                                'Error  Longitud(Km)':'ERROR LONGITUD (Km)',
                                'Error  Profundidad(Km)':'ERROR PROFUNDIDAD (Km)'}, inplace = True)
df_sismos.drop(['Fecha  (UTC)','Hora  (UTC)'],axis=1,inplace=True)
df_sismos=df_sismos[(df_sismos['PROF. (Km)']<=32)&(df_sismos['PROF. (Km)']>(z_topo.min()*(-1/1000)))&(df_sismos['MAGNITUD']>0)& #Filtros previos
        (df_sismos['LATITUD (°)']>lai)&(df_sismos['LATITUD (°)']<las)
        &(df_sismos['LONGITUD (°)']>loi)&(df_sismos['LONGITUD (°)']<los)]
df_sismos['PROF. (m)']=-df_sismos['PROF. (Km)']*1000 #
df_sismos['ERROR PROFUNDIDAD (m)']=df_sismos['ERROR PROFUNDIDAD (Km)']*1000 #Conversion de km a m
df_sismos['ERROR LONGITUD (°)']=df_sismos['ERROR LONGITUD (Km)']/111.1 #Conversion de km a °
df_sismos['ERROR LATITUD (°)']=df_sismos['ERROR LATITUD (Km)']/111.1 
df_sismos['FECHA - HORA UTC']=df_sismos['FECHA - HORA UTC'].apply(lambda x:pd.to_datetime(x))#Conversion de str a UTCDateTime

#Errores con topografia
df_sismos_err=df_sismos[df_sismos['PROF. (m)']+df_sismos['ERROR PROFUNDIDAD (m)']>(z_topo.min())]
df_sismos_no_err=df_sismos[df_sismos['PROF. (m)']+df_sismos['ERROR PROFUNDIDAD (m)']<=(z_topo.min())]
df_sismos_no_err.loc[:,'ERROR PROFUNDIDAD SUP (m)']=df_sismos_no_err.loc[:,'ERROR PROFUNDIDAD (m)'].copy()
z=[]
for x,y,zhip,zerr in zip(df_sismos_err['LONGITUD (°)'],df_sismos_err['LATITUD (°)'],
                df_sismos_err['PROF. (m)'],df_sismos_err['ERROR PROFUNDIDAD (m)']):
    df_elev=df_topo[(df_topo[0]<(x+0.005))&(df_topo[0]>(x-0.005))&
                (df_topo[1]<(y+0.005))&(df_topo[1]>(y-0.005))]
    d=1
    for x0,y0,z0 in zip(df_elev[0],df_elev[1],df_elev[2]):
        dist=np.sqrt(((x-x0)**2)+((y-y0)**2))
        if dist<d:
            d=dist
            z1=z0
    if z1<=(zhip+zerr):
        z.append(z1-zhip)
    else :
        z.append(zerr)
df_sismos_err.loc[:,'ERROR PROFUNDIDAD SUP (m)']=(np.array(z))
df_sismos=pd.concat([df_sismos_err, df_sismos_no_err])
del df_sismos_err
del df_sismos_no_err

#Inyeccion de H2O
iny=pd.read_csv(r'datasets\inyeccion_geo.csv',delimiter=';',decimal=',')
iny=iny[:-1]

inyec=[]
for name,lon,lat,alt in zip(iny['CAMPO'].apply(lambda x:str(x)),iny['X'],iny['Y'],[2100]*len(iny['CAMPO'])):
    un=dict(
            showarrow=False,
            x=lon,
            y=lat,
            z=alt+10,
            text=name,
            xanchor="left",
            xshift=10,
            opacity=0.7,
            font=dict(
                color="black",
                size=12
            ))
    inyec.append(un)

#Kale
df_kale=pd.read_csv('datasets/kale.csv')
df_kale['msnm']=[69]*3

kale= go.Scatter3d(
    x=np.array(-73.85660000000),
    y=np.array(7.36551000000),
    z=np.array(69+100),
    mode='markers',
    marker_symbol='diamond',
    name="PPII Kalé",
    hovertemplate ='PPII Kalé',
    marker=dict(
        size=10,
        color='gold'
    )
)
#Semaforo sismico Kale
pozo_inv_kale = df_kale[df_kale['Tipo']=='Investigación']
x_pozo_inv_kale, y_pozo_inv_kale  = pozo_inv_kale['Longitud'].values[0], pozo_inv_kale['Latitud'].values[0]
h_pozo_inv_kale = 3.902 #km
h_pozo_inv_kale_m = h_pozo_inv_kale*1000 #m
r_ext = 2*h_pozo_inv_kale_m+20000 #m

#Asignamos las dimensiones y ubicacion del cilindro interno y externo respectivamente
r1 = 2*h_pozo_inv_kale_m /(111.1*1000) #Radio interno es dos veces la profundidad medida del pozo. De acuerdo con Resolución 40185 del 2020 del MME. Profundidad aproximada en pozo de investigación es 3902 m
a1 = 0 #Altura
h1 = -16000 #Profundidad del cilindro de 16 km
x01=float(x_pozo_inv_kale)
y01=float(y_pozo_inv_kale)

r2 = (2*h_pozo_inv_kale_m+20000)/(111.1*1000) #Radio externo es  2*h (profundidad del pozo) + 20 km
a2 = 0 #Altura
h2 = -16000 #Profundidad del cilindro de 16 km
x02=float(x_pozo_inv_kale)
y02=float(y_pozo_inv_kale)

r3 = (50000)/(111.1*1000) #Radio externo es  2*h (profundidad del pozo) + 20 km
a3 = 0 #Altura
h3 = -32000 #Profundidad del cilindro de 16 km
x03=float(x_pozo_inv_kale)
y03=float(y_pozo_inv_kale)

#Efectuamos los calculos correspondientes  a la funcion
x1, y1, z1 = cylinder(r1, h1,x01,y01, a=a1)
x2, y2, z2 = cylinder(r2, h2,x02,y02, a=a2)
x3, y3, z3 = cylinder(r3, h3,x03,y03, a=a3)

#Elaboramos la proyeccion para el volumen de suspension
cyl1 = go.Surface(x=x1, y=y1, z=z1,
                 colorscale = [[0, 'red'],[1, 'red']],#El color se da porque aqui es donde se analizan los sismos 
                                                            #que pueden dar un alarma verde,amarilla o naranja
                 showscale=False,
                 opacity=0.5,
                 name='Volumen monitoreo estado rojo')
xb_low, yb_low, zb_low = boundary_circle(r1, a1,x01,y01)
xb_up, yb_up, zb_up = boundary_circle(r1, a1+h1,x01,y01)

bcircles1 =go.Scatter3d(x = xb_low.tolist()+[None]+xb_up.tolist(),
                        y = yb_low.tolist()+[None]+yb_up.tolist(),
                        z = zb_low.tolist()+[None]+zb_up.tolist(),
                        mode ='lines',
                        line = dict(color='red', width=2),
                        opacity =0.55, showlegend=False,
                        name='Volumen monitoreo estado rojo')

#Elaboramos la proyeccion para el volumen de monitoreo
cyl2 = go.Surface(x=x2, y=y2, z=z2,
                 colorscale = [[0, 'green'],[1, 'orange']],
                 showscale=False,
                 opacity=0.7,
                 name='Volumen monitoreo para estado verde,amarillo y naranja')

xb_low, yb_low, zb_low = boundary_circle(r2, a2,x02,y02)
xb_up, yb_up, zb_up = boundary_circle(r2,a2+h2,x02,y02)

#Bordes
bcircles2 =go.Scatter3d(x = xb_low.tolist()+[None]+xb_up.tolist(),
                        y = yb_low.tolist()+[None]+yb_up.tolist(),
                        z = zb_low.tolist()+[None]+zb_up.tolist(),
                        mode ='lines',
                        line = dict(color='green', width=2),
                        opacity =0.75, showlegend=False,
                        name='Volumen monitoreo para estado verde,amarillo y naranja'
                        )

#Elaboramos la proyeccion para el cilindro de volumen externo
cyl3 = go.Surface(x=x3, y=y3, z=z3,
                 colorscale = [[0, 'aqua'],[1, 'aqua']],
                 showscale=False,
                 opacity=0.4,
                 name='Volumen externo')

xb_low, yb_low, zb_low = boundary_circle(r3, a3,x03,y03)
xb_up, yb_up, zb_up = boundary_circle(r3,a3+h3,x03,y03)

#Bordes
bcircles3 =go.Scatter3d(x = xb_low.tolist()+[None]+xb_up.tolist(),
                        y = yb_low.tolist()+[None]+yb_up.tolist(),
                        z = zb_low.tolist()+[None]+zb_up.tolist(),
                        mode ='lines',
                        line = dict(color='blue', width=2),
                        opacity =0.75, showlegend=False,
                        name='Volumen externo'
                        )

#Platero
platero= go.Scatter3d(
    x=np.array(-73.8938980),
    y=np.array(7.2572498),
    z=np.array(69+100),
    mode='markers',
    marker_symbol='diamond',
    name="PPII Platero",
    hovertemplate ='PPII Platero',
    marker=dict(
        size=10,
        color='gold'
    )
)
#Semaforo sismico Platero
x_pozo_inv_plat, y_pozo_inv_plat  = -73.8938980,7.2572498
h_pozo_inv_plat_m = 3227.8 #m
r_ext_plat = 2*h_pozo_inv_plat_m+20000 #m

#Asignamos las dimensiones y ubicacion del cilindro interno y externo respectivamente
r1 = 2*h_pozo_inv_plat_m /(111.1*1000) #Radio interno es dos veces la profundidad medida del pozo. De acuerdo con Resolución 40185 del 2020 del MME. Profundidad aproximada en pozo de investigación es 3902 m
a1 = 0 #Altura
h1 = -16000 #Profundidad del cilindro de 16 km
x01=float(x_pozo_inv_plat)
y01=float(y_pozo_inv_plat)

r2 = (2*h_pozo_inv_plat_m+20000)/(111.1*1000) #Radio externo es  2*h (profundidad del pozo) + 20 km
a2 = 0 #Altura
h2 = -16000 #Profundidad del cilindro de 16 km
x02=float(x_pozo_inv_plat)
y02=float(y_pozo_inv_plat)

r3 = (50000)/(111.1*1000) #Radio externo es  2*h (profundidad del pozo) + 20 km
a3 = 0 #Altura
h3 = -32000 #Profundidad del cilindro de 16 km
x03=float(x_pozo_inv_plat)
y03=float(y_pozo_inv_plat)

#Efectuamos los calculos correspondientes  a la funcion
x1, y1, z1 = cylinder(r1, h1,x01,y01, a=a1)
x2, y2, z2 = cylinder(r2, h2,x02,y02, a=a2)
x3, y3, z3 = cylinder(r3, h3,x03,y03, a=a3)

#Elaboramos la proyeccion para el volumen de suspension
cyl1p = go.Surface(x=x1, y=y1, z=z1,
                 colorscale = [[0, 'red'],[1, 'red']],#El color se da porque aqui es donde se analizan los sismos 
                                                            #que pueden dar un alarma verde,amarilla o naranja
                 showscale=False,
                 opacity=0.5,
                 name='Volumen monitoreo estado rojo')
xb_low, yb_low, zb_low = boundary_circle(r1, a1,x01,y01)
xb_up, yb_up, zb_up = boundary_circle(r1, a1+h1,x01,y01)

bcircles1p =go.Scatter3d(x = xb_low.tolist()+[None]+xb_up.tolist(),
                        y = yb_low.tolist()+[None]+yb_up.tolist(),
                        z = zb_low.tolist()+[None]+zb_up.tolist(),
                        mode ='lines',
                        line = dict(color='red', width=2),
                        opacity =0.55, showlegend=False,
                        name='Volumen monitoreo estado rojo')

#Elaboramos la proyeccion para el volumen de monitoreo
cyl2p = go.Surface(x=x2, y=y2, z=z2,
                 colorscale = [[0, 'green'],[1, 'orange']],
                 showscale=False,
                 opacity=0.7,
                 name='Volumen monitoreo para estado verde,amarillo y naranja')

xb_low, yb_low, zb_low = boundary_circle(r2, a2,x02,y02)
xb_up, yb_up, zb_up = boundary_circle(r2,a2+h2,x02,y02)

#Bordes
bcircles2p =go.Scatter3d(x = xb_low.tolist()+[None]+xb_up.tolist(),
                        y = yb_low.tolist()+[None]+yb_up.tolist(),
                        z = zb_low.tolist()+[None]+zb_up.tolist(),
                        mode ='lines',
                        line = dict(color='green', width=2),
                        opacity =0.75, showlegend=False,
                        name='Volumen monitoreo para estado verde,amarillo y naranja'
                        )

#Elaboramos la proyeccion para el cilindro de volumen externo
cyl3p = go.Surface(x=x3, y=y3, z=z3,
                 colorscale = [[0, 'aqua'],[1, 'aqua']],
                 showscale=False,
                 opacity=0.4,
                 name='Volumen externo')

xb_low, yb_low, zb_low = boundary_circle(r3, a3,x03,y03)
xb_up, yb_up, zb_up = boundary_circle(r3,a3+h3,x03,y03)

#Bordes
bcircles3p =go.Scatter3d(x = xb_low.tolist()+[None]+xb_up.tolist(),
                        y = yb_low.tolist()+[None]+yb_up.tolist(),
                        z = zb_low.tolist()+[None]+zb_up.tolist(),
                        mode ='lines',
                        line = dict(color='blue', width=2),
                        opacity =0.75, showlegend=False,
                        name='Volumen externo'
                        )

#Estaciones sismologicas
df_sta_vmm=pd.read_csv('datasets/VMM_STA.csv',delimiter=';',decimal=',')
df_sta_lom=pd.read_csv('datasets//LOMA_STA.csv',delimiter=';',decimal=',')

STA_VMM = go.Scatter3d(
    x=df_sta_vmm['LONGITUD'],
    y=df_sta_vmm['LATITUD'],
    z=df_sta_vmm['ALTITUD (msnm)']+10, #Para sobresalir de la topografía
    mode='markers',
    marker_symbol='diamond',
    name="Estación sismologica VMM",
    hovertemplate ='Longitud:'+df_sta_vmm['LONGITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Latitud:'+df_sta_vmm['LATITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Elevacion:'+df_sta_vmm['ALTITUD (msnm)'].apply(lambda x:str(x))+'msnm <br>'+
                    'Nombre de la estación:'+df_sta_vmm['NOMBRE ESTACIÓN']+'°'+'<br>'+
                    'Código:'+df_sta_vmm['CODIGO']+'<br>'+
                    'Agencia:'+df_sta_vmm['AGENCIA']+'<br>'+
                    'Fecha de instalación:'+df_sta_vmm['FECHA DE INSTALACIÓN'].apply(lambda x:str(x))+'<br>'+
                    'Fecha de retiro:'+df_sta_vmm['FECHA DE RETIRO'].apply(lambda x:str(x))+'<br>'+
                    'Estado:'+df_sta_vmm['ESTADO'],
    marker=dict(
        size=4,
        color='blueviolet'
    )
)
STA_LOM = go.Scatter3d(
    x=df_sta_lom['LONGITUD'],
    y=df_sta_lom['LATITUD'],
    z=df_sta_lom['ALTITUD (msnm)']+10,
    mode='markers',
    marker_symbol='diamond',
    name="Estación sismologica la Loma, Cesar",
    hovertemplate ='Longitud:'+df_sta_lom['LONGITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Latitud:'+df_sta_lom['LATITUD'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Elevacion:'+df_sta_lom['ALTITUD (msnm)'].apply(lambda x:str(x))+'msnm <br>'+
                    'Nombre de la estación:'+df_sta_lom['NOMBRE ESTACIÓN']+'°'+'<br>'+
                    'Código:'+df_sta_lom['CODIGO']+'<br>'+
                    'Agencia:'+df_sta_lom['AGENCIA']+'<br>'+
                    'Fecha de instalación:'+df_sta_lom['FECHA DE INSTALACIÓN'].apply(lambda x:str(x))+'<br>'+
                    'Fecha de retiro:'+df_sta_lom['FECHA DE RETIRO'].apply(lambda x:str(x))+'<br>'+
                    'Estado:'+df_sta_lom['ESTADO'],
    marker=dict(
        size=6,
        color='blueviolet'
    )
)

#Rios
df_rivers=pd.read_csv('datasets\drenajes_dem15arcs_WGS84_LBG.csv')

#Cargar datos de pozos
df_pozos=pd.read_csv('datasets/pozos.csv',usecols=['lon', 'lat', 'UWI', 'WELL_NAME', 
'DEPARTAMEN', 'WELL_COU_1', 'WELL_TVD', 'WELL_KB_EL',
       'ROTARY_ELE', 'WELL_DRILL', 'WELL_GROUN', 'FIELD_ABRE',
       'CONTRATO', 'WELL_SPUD_', 'COORD_QUAL', 'COMMENT_', 'WELL_COMPL', 'WELL_STA_1',
       'WELLTYPE', 'FECHA_ACTU',
       'OPERATOR_W', 'COMPANY_CO', 'z'])
Pozos = go.Scatter3d(
    x=df_pozos['lon'],
    y=df_pozos['lat'],
    z=df_pozos['z']+15,
    mode='markers',
    name="Pozo petrolífero",
    hovertemplate ='Longitud:'+df_pozos['lon'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Latitud:'+df_pozos['lat'].apply(lambda x:str(x))+'°'+'<br>'+
                    'Elevacion:'+df_pozos['z'].apply(lambda x:str(x))+'msnm <br>'+
                    'UWI:'+df_pozos['UWI']+'°'+'<br>'+
                    'Nombre del pozo:'+df_pozos['WELL_NAME']+'<br>'+
                    'Departamento:'+df_pozos['DEPARTAMEN']+'<br>'+
                    'Municipio:'+df_pozos['WELL_COU_1']+'<br>'+
                    'Tipo:'+df_pozos['WELLTYPE']+'<br>'+
                    'Operador:'+df_pozos['OPERATOR_W']+'<br>'+
                    'Compañia:'+df_pozos['COMPANY_CO'],
    marker=dict(
        size=1.5,
        color='black'
    )
)

#Rezumaderos
df_rezumaderos=pd.read_csv('datasets\REZUMADEROS_WGS84_LBG.txt',decimal=',',delimiter=';')
rez_txt=('Longitud:'+df_rezumaderos['X'].apply(lambda x:str(x))+'°'+
            '<br>'+'Latitud:'+df_rezumaderos['Y'].apply(lambda x:str(x))+'°'+
            '<br>'+'Elevacion:'+df_rezumaderos['Z'].apply(lambda x:str(x))+'msnm'+
            '<br>'+'Tipo:'+df_rezumaderos['TIPO'].apply(lambda x:str(x))+
            '<br>'+'Autor:'+df_rezumaderos['AUTOR'].apply(lambda x:str(x))+
            '<br>'+'Empresa:'+df_rezumaderos['EMPRESA'].apply(lambda x:str(x))+
            '<br>'+'Formación:'+df_rezumaderos['FORMACION_'].apply(lambda x:str(x))+
            '<br>'+'Tipo secundario:'+df_rezumaderos['TIPO_2'].apply(lambda x:str(x))+
            '<br>'+'Departamento:'+df_rezumaderos['DPTOS_NOMB'].apply(lambda x:str(x))+
            '<br>'+'Capital:'+df_rezumaderos['CAPITAL'].apply(lambda x:str(x)))
rez = go.Scatter3d(
    x=df_rezumaderos['X'],
    y=df_rezumaderos['Y'],
    z=df_rezumaderos['Z'],
    mode='markers',
    name="Rezumaderos",
    hovertemplate =rez_txt,
    marker=dict(
        size=5,
        color='gray'
    )
)
#Poblaciones
df_poblaciones=pd.read_csv('datasets/poblaciones.csv',usecols=['Name','Longitud','Latitud','SAMPLE_1'])
Poblaciones = go.Scatter3d(
    x=df_poblaciones['Longitud'],
    y=df_poblaciones['Latitud'],
    z=df_poblaciones['SAMPLE_1']+10, #Conseguir alturas
    mode='markers',
    name="Población",
    marker_symbol='square',
    hovertemplate =df_poblaciones['Name'],
    marker=dict(
        size=6,
        color='red'
    ),
    textposition="bottom right"
)
Pobl=[]
for name,lon,lat,alt in zip(df_poblaciones['Name'],df_poblaciones['Longitud'],df_poblaciones['Latitud'],df_poblaciones['SAMPLE_1']):
    un=dict(
            showarrow=False,
            x=lon,
            y=lat,
            z=alt+10,
            text=name,
            xanchor="left",
            xshift=10,
            opacity=0.7,
            font=dict(
                color="black",
                size=12
            ))
    Pobl.append(un)

#Carreteras
viass=pd.read_csv('datasets\Vias_LBG.csv')

#Fallas
faults=pd.read_csv('datasets\Fallas_LBG.csv')
faults=faults.drop('Unnamed: 0',axis=1)

#Campos
campet=pd.read_csv('datasets\campos_LBG.csv',decimal=',')
campet['X']=campet['X'].apply(lambda x:float(x))
campet['Y']=campet['Y'].apply(lambda x:float(x))
campet['Z']=campet['Z'].apply(lambda x:float(x))
campet_1=pd.read_csv('datasets/campos_1_LBG.csv')
campet_1=campet_1.drop_duplicates(subset=['LINE_ID'])

# ls_x_s,ls_y_s,ls_z_s=lin_list('datasets/lineas.csv') #Lineas sismicas
linsis=pd.read_csv('datasets\lineas_LBG.csv',decimal=',')
linsis['X']=linsis['X'].apply(lambda x:float(x))
linsis['Y']=linsis['Y'].apply(lambda x:float(x))
linsis['Z']=linsis['Z'].apply(lambda x:float(x))
linsis_1=pd.read_csv('datasets/lineas_1_LBG.csv')
linsis_1=linsis_1.drop_duplicates(subset=['LINE_ID'])



Eoceno=geology('datasets/DISCORDANCIA_EOCENO.txt','pink','red','Discordancia del Eoceno Medio')
Colorado=geology('datasets/TOPE_COLORADO.txt','darkblue','aquamarine','Tope Formación Colorado')
Mugrosa=geology('datasets/TOPE_MUGROSA.txt','green','greenyellow','Tope Formación Mugrosa')
Chorros=geology('datasets/TOPE_CHORROS.txt','orangered','yellow','Tope Grupo Chorros')
Real=geology('datasets/BASE_CUATERNARIO.txt','purple','pink','Tope Grupo Real')





SISMICA=img_3d("assets/perfil_2.jpg",-74.115,7.58,-72.954,6.806,4300,-20000)


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SUPERHERO])
app.config['suppress_callback_exceptions'] = True
#Cargars los datos

card_main=dbc.Card(
    dbc.CardBody(
        
            [dbc.Nav([
                dbc.NavLink("Inicio", href="https://www.centrodetransparenciappii.org/", active="exact"),
                dbc.NavLink("Modelo 3D VMM", href="", active="exact"),
                dbc.NavLink("Semáforo sísmico", href="https://pinguinodigital.com/wp-content/uploads/2020/08/pagina-en-construcci%C3%B3n1.jpg", active="exact"),
            ]),
            html.H2("Modelo Tridimensional de Sismicidad en el Valle Medio del Magdalena", className="card-title"),
            html.H4("Transparencia de la superficie:", className="card-subtitle"),
            dcc.Slider(
                id='TOPO',
                min=0,
                max=1,
                step=0.1,
                value=1,
                tooltip={"placement": "bottom", "always_visible": True}),
            html.H4("Exageración vertical:", className="card-subtitle"),
            dcc.Slider(
                id='EXG',
                min=10,
                max=100,
                step=10,
                value=30,
                tooltip={"placement": "bottom", "always_visible": True}),
            html.H4("Magnitudes:", className="card-subtitle"),
                    dcc.RangeSlider(
                id='MAGN',
                min=df_sismos['MAGNITUD'].min(),
                max=df_sismos['MAGNITUD'].max(),
                step=0.1,
                value=[3, df_sismos['MAGNITUD'].max()],
                allowCross=False,
                tooltip={"placement": "bottom", "always_visible": True}
            ),
            html.H4("Profundidad (m):", className="card-subtitle"),
            dcc.RangeSlider(
                id='DEPTH',
                min=df_sismos['PROF. (m)'].min(),
                max=df_sismos['PROF. (m)'].max(),
                step=100,
                value=[df_sismos['PROF. (m)'].min(), df_sismos['PROF. (m)'].max()],
                allowCross=False,
                tooltip={"placement": "bottom", "always_visible": True}
            ),
            html.H4("Fecha:", className="card-subtitle"),
            dcc.DatePickerRange(
                    id='DATE',
                    start_date_placeholder_text="Start Date",
                    end_date_placeholder_text="End Date",
                    calendar_orientation='horizontal',
                    start_date=df_sismos['FECHA - HORA UTC'].min(),
                    end_date=df_sismos['FECHA - HORA UTC'].max(),
                    day_size=30,
                    min_date_allowed=df_sismos['FECHA - HORA UTC'].min(),
                    max_date_allowed=df_sismos['FECHA - HORA UTC'].max(),
                    persistence=True,
                    #initial_visible_month=df_sismos['FECHA - HORA UTC'].min(),
                    reopen_calendar_on_clear=False
                ),
                # html.H4("Perfiles:", className="card-subtitle"),
                # html.H5("Punto 1 (Longitud-Latitud)", className="card-subtitle"),
                # dcc.Input(id="Longitud 1", type="number", placeholder="Longitud 1", min=loi, max=los, step=0.01,style={'marginRight':'10px'},value=loi),
                # dcc.Input(id="Latitud 1", type="number", placeholder="Latitud 1", min=lai, max=las, step=0.01, debounce=True,value=lai),
                # html.H5("Punto 2 (Longitud-Latitud)", className="card-subtitle"),
                # dcc.Input(id="Longitud 2", type="number", placeholder="Longitud 2", min=loi, max=los, step=0.01,value=los,style={'marginRight':'10px'}),
                # dcc.Input(id="Latitud 2", type="number", placeholder="Latitud 2", min=lai, max=las, step=0.01,value=las,debounce=True),
            html.H4("Variables de sismicidad:", className="card-subtitle"),
            dcc.Dropdown(id='SEISMO',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': 'Localización', 'value': 'LOC'},
                            {'label': 'Fecha', 'value': 'FEC'},
                            {'label': 'Magnitud', 'value': 'MAG'},
                            {'label': 'RMS', 'value': 'RMS'},
                            {'label': 'Errores', 'value': 'ERR'},
                            
                        ],
                        value=['LOC', 'FEC','MAG','RMS','ERR'],
                        multi=True
                    ),
            html.H4("________________________________________", className="card-subtitle"),
            html.H4("Cartografía y linea base:", className="card-subtitle"),
            dcc.Dropdown(id='CART',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': ' Pozo Kalé (ANH)', 'value': 'KALE'},
                            {'label': ' Cilindro en suspensión Semáforo sísmico para Kalé (SGC)', 'value': 'SEM_KALE'},
                            {'label': ' Pozo Platero (ANH)', 'value': 'PLATERO'},
                            {'label': ' Cilindro en suspensión Semáforo sísmico para Platero (SGC)', 'value': 'SEM_PLATERO'},
                            {'label': ' Barras de Error (SGC)', 'value': 'ERROR'},
                            {'label': ' Estaciones sismológicas (SGC)', 'value': 'STA'},
                            {'label': ' Poblaciones (UNAL-ANH-MINCIENCIAS)', 'value': 'POB'},
                            {'label': ' Drenajes (IGAC)', 'value': 'RIV'},
                            {'label': ' Vias (IGAC)', 'value': 'VIA'},
                            
                        ],
                        value=[],
                        multi=True
                    ),
            html.H4("________________________________________", className="card-subtitle"),
            html.H4("Información petrolífera y sísmica:", className="card-subtitle"),
            dcc.Dropdown(id='PETRO',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': ' Pozos petrolíferos (UNAL-ANH-MINCIENCIAS)', 'value': 'POZO'},
                            {'label': ' Campos petrolíferos (UNAL-ANH-MINCIENCIAS)', 'value': 'FIELD'},
                            {'label': ' Trazo en superficie de líneas sísmicas (UNAL-ANH-MINCIENCIAS)', 'value': 'LIN'},
                            {'label': ' Rezumaderos (ANH)', 'value': 'REZ'},
                            {'label': ' ANH-TR-2006-04-A (ANH)', 'value': 'SEIS'}

                            
                        ],
                        value=[],
                        multi=True
                    ),
            html.H4("Inyección de Agua:", className="card-subtitle"),
            dcc.Dropdown(id='INY',
                                    placeholder="Fecha inyeccion",
                                    style={'color': 'black'},
                                    options=[{'label':x,'value':x} for x in iny.columns[6:]],
                                    value=[],
                                    multi=False
                                ),
            html.H4("________________________________________", className="card-subtitle"),
            html.H4("Geología:", className="card-subtitle"),
            dcc.Dropdown(id='GEOL',
                        placeholder="Variables a desplegar...",
                        style={'color': 'black'},
                        options=[
                            {'label': ' Fallas Geológicas (SGC)', 'value': 'FALL'},
                            {'label': ' Tope Grupo Real (UNAL-ANH-MINCIENCIAS)', 'value': 'REAL'},
                            {'label': ' Tope Formación Colorado (UNAL-ANH-MINCIENCIAS)', 'value': 'COL'},
                            {'label': ' Tope Formación Mugrosa (UNAL-ANH-MINCIENCIAS)', 'value': 'MUG'},
                            {'label': ' Tope Grupo Chorros (UNAL-ANH-MINCIENCIAS)', 'value': 'CHO'},
                            {'label': ' Discordancia del Eoceno Medio (UNAL-ANH-MINCIENCIAS)', 'value': 'EOC'},
                            # {'label': ' Geología superficial (SGC)', 'value': 'GEO'},
                            
                        ],
                        value=[],
                        multi=True
                    ),
            dcc.Markdown('''
                    * ANH: Agencia Nacional de Hidrocarburos
                    * MINCIENCIAS: Ministerio de Ciencia Tecnología e Innovación
                    * SGC: Servicio Geológico Colombiano
                    * UNAL: Universidad Nacional de Colombia
                '''),html.A(
                                html.Button("Descarga como HTML"), 
                                id="download",
                                download="VMM_SEISMO.html"
                            ),
                     dbc.CardImg(src="assets\logos.png", bottom=True, alt='Logos_convenio_tripartito',)    
                 ,], 
        ),
    color="secondary",   # https://bootswatch.com/default/ for more card colors
    inverse=True,   # change color of text (black or white)
    # outline=False,  # True = remove the block colors from the background and header,
    
)


card_references=dbc.Card(
    dbc.CardBody([
        html.H2("Referencias", className="card-title"),
        html.H6("Agencia Nacional de Hidrocarburos - ANH & Servicio Geológico Colombiano - SGC (2016). Informe final del Convenio interadministrativo 194 ANH-014 SGC, entre la Agencia Nacional de Hidrocarburos y el Servicio Geológico Colombiano.", 
            className="card-text"),
        html.H6("Agencia Nacional de Hidrocarburos - ANH (2010). Mapa de Rezumaderos. Información Geológica y Geofísica. https://www.anh.gov.co/Informacion-Geologica-y-Geofisica/Estudios-Integrados-y-Modelamientos/Paginas/MAPA-DE-REZUMADEROS.aspx", 
            className="card-text"),         
        html.H6("Ángel-Martínez, C.E., Prieto-Gómez, G.A., Cristancho-Mejía, F., Sarmiento-Orjuela, A.M., Vargas-Quintero, J.A., Delgado-Mateus, C.J., Torres-Rojas, E., Castelblanco-Ossa, C.A., Camargo-Rache, G.L., Amazo-Gómez, D.F., Cipagauta-Mora, J.B., Lucuara-Reyes, E.D., Ávila-López, K.L. Fracica-González, L.R., Martín-Ravelo, A.S., Atuesta-Ortiz, D.A., Gracía-Romero, D.F., Triviño Cediel , R.J., Jaimes Villarreal, V.N., y Alarcón Rodríguez, W.F.(2021). Proyecto MEGIA: Modelo Geológico-Geofísico del Valle Medio del Magdalena. Producto No. 5. Bogotá: 192 pp.", 
            className="card-text"),
        html.H6("Dionicio, V., Mercado, O. y Lizarazo, M. (2020). Semáforo para el monitoreo sísmico durante el desarrollo de los proyectos piloto de investigación integral en yacimientos no convencionales de hidrocarburos en Colombia. Bogotá: Servicio Geológico Colombiano.", 
            className="card-text"),
        html.H6("Gómez, J. & Montes, N.E., compiladores. 2020. Mapa Geológico de Colombia 2020. Escala 1:1 000 000. Servicio Geológico Colombiano, 2 hojas. Bogotá.​", 
            className="card-text"),
        html.H6("Instituto Geográfico Agustin Codazzi - IGAC (2019). Base de datos vectorial básica. Colombia. Escala 1:100.000. Colombia en Mapas. https://www.colombiaenmapas.gov.co/#", 
            className="card-text"),
        html.H6("Servicio Geológico Colombiano. (2021). Banco de Información Petrolera. https://srvags.sgc.gov.co/JSViewer/GEOVISOR_BIP/", 
            className="card-text"),
        html.H6("Servicio Geológico Colombiano. (2022). Catálogo línea base de sismicidad: Valle Medio del Magdalena y La Loma Cesar. http://bdrsnc.sgc.gov.co/paginas1/catalogo/Consulta_Valle_Medio/valle_medio.php", 
            className="card-text"),
         html.H6("Tozer, B, Sandwell, D. T., Smith, W. H. F., Olson, C., Beale, J. R., & Wessel, P. (2019). Global bathymetry and topography at 15 arc sec: SRTM15+. Distributed by OpenTopography. https://doi.org/10.5069/G92R3PT9. Accessed: 2022-02-10", 
            className="card-text"),
    ]))

app.layout = html.Div([
    dbc.Row([dbc.Col(card_main, width=12),
             dbc.Col(card_references, width=12)], 
             justify="start"), 
             
])


@app.callback(
     [dash.dependencies.Output(component_id='download', component_property='href'),
      dash.dependencies.Output(component_id='DATE', component_property='initial_visible_month')],



    [dash.dependencies.Input(component_id='TOPO', component_property='value'),
     dash.dependencies.Input(component_id='EXG', component_property='value'),
     dash.dependencies.Input(component_id='DATE', component_property='start_date'),
     dash.dependencies.Input(component_id='DATE', component_property='end_date'),
     dash.dependencies.Input(component_id='MAGN', component_property='value'),
     dash.dependencies.Input(component_id='DEPTH', component_property='value'),
     dash.dependencies.Input(component_id='SEISMO', component_property='value'),
     dash.dependencies.Input(component_id='CART', component_property='value'),
     dash.dependencies.Input(component_id='PETRO', component_property='value'),
     dash.dependencies.Input(component_id='INY', component_property='value'),
     dash.dependencies.Input(component_id='GEOL', component_property='value'),
    #  dash.dependencies.Input(component_id='Longitud 1', component_property='value'),
    #  dash.dependencies.Input(component_id='Longitud 2', component_property='value'),
    #  dash.dependencies.Input(component_id='Latitud 1', component_property='value'),
    #  dash.dependencies.Input(component_id='Latitud 2', component_property='value') 
     ])

def update_figure(TOPO,EXG,START_DATE,END_DATE,MAGN,DEPTH,SEISMO,CART,PETRO,INY,GEOL):
        fig=go.Figure()
        # if np.isin('GEO', GEOL):
        #     if TOPO==0:
        #         DISM=0
        #     else:
        #         DISM=0.01
        #     fig.add_trace(go.Surface(z=z_topo,showscale=False, x=x_topo, y=y_topo,colorscale=['black','black'],lighting=dict(ambient=0.3,diffuse=0.5),
        #             showlegend=False,opacity=TOPO-DISM,name='Topografía'))
        #     directory = 'datasets\CSV_UNIDADES'
        #     for filename in os.scandir(directory):
        #         if filename.is_file():
        #             name=(str(filename.path).split('\\'))[-1]
        #             name=name.replace('.csv','')
        #             name=name.replace('_','?')
        #             df_1=df_new[df_new['SimboloUC']==name]
        #             text='Edad: '+np.array(df_1['Edad'].apply(lambda x:str(x)))[0]+'<br>Descripción: '+np.array(df_1['Descripcio'].apply(lambda x:str(x)))[0]#+'<br>UGIntegrad: '+np.array(df_1['UGIntegrad'].apply(lambda x:str(x)))[0]
        #             text=str(text)
        #             fig.add_trace(geology_super(filename.path,np.array(df_1['Color'].apply(lambda x:str(x)))[0],name,text,TOPO))
            
        # else:
        fig.add_trace(go.Surface(z=z_topo,showscale=False, x=x_topo, y=y_topo,colorscale=['green','greenyellow','saddlebrown','saddlebrown','saddlebrown','saddlebrown','snow','snow'],
                    showlegend=False,opacity=TOPO,name='Topografía',lighting=dict(ambient=0.3,diffuse=0.5)))
        df_sismos_1=df_sismos[(df_sismos['FECHA - HORA UTC']<=END_DATE)&(df_sismos['FECHA - HORA UTC']>=START_DATE)&
        (df_sismos['MAGNITUD']>=MAGN[0])&(df_sismos['MAGNITUD']<=MAGN[1])
        &(df_sismos['PROF. (m)']>=DEPTH[0])&(df_sismos['PROF. (m)']<=DEPTH[1])]
        text=text_scatter(SEISMO,df_sismos_1)
        if np.isin('ERROR', CART):
            err=True
        else:
            err=False
        fig.add_trace(go.Scatter3d(
            x=df_sismos_1['LONGITUD (°)'],y=df_sismos_1['LATITUD (°)'],z=df_sismos_1['PROF. (m)'],mode='markers',
            marker=dict(
                size=(df_sismos_1['MAGNITUD'])*3,
                color=df_sismos_1['PROF. (m)'],                # set color to an array/list of desired values
                colorscale='Jet',   # choose a colorscale
                opacity=0.8,
                cmax=df_sismos['PROF. (m)'].max(),
                cmin=-32000,
            ),
            error_x=dict(
                array=df_sismos_1['ERROR LONGITUD (°)'],                # set color to an array/list of desired values
                color='red',   # choose a colorscale
                symmetric=True,
                width=0.01,
                visible=err
            ),
            error_y=dict(
                array=df_sismos_1['ERROR LATITUD (°)'],                # set color to an array/list of desired values
                color='red',   # choose a colorscale
                symmetric=True,
                width=0.01,
                visible=err
            ),
            error_z=dict(
                array=df_sismos_1['ERROR PROFUNDIDAD SUP (m)'], 
                arrayminus=df_sismos_1['ERROR PROFUNDIDAD (m)'] ,             
                color='red',   # choose a colorscale
                symmetric=False,
                width=0.01,
                visible=err
            ),
            hovertemplate=text,
                name='Sismos',
                showlegend=False))
        if np.isin('KALE', CART):
            fig.add_trace(kale)
        if np.isin('RIV', CART):
            for i in df_rivers['DRENAJE'].unique():
                riv=df_rivers[df_rivers['DRENAJE']==i]
                fig.add_trace(go.Scatter3d(z=riv['Z'], x=riv['X'], y=riv['Y'],mode='markers',
                name=str(i),marker_symbol='square',marker=dict(color='aqua',size=3)))
        if np.isin('SEM_KALE', CART):
            fig.add_traces(data=[cyl1, bcircles1,cyl2, bcircles2,cyl3,bcircles3])
        if np.isin('PLATERO', CART):
            fig.add_trace(platero)
        if np.isin('SEM_PLATERO', CART):
            fig.add_traces(data=[cyl1p, bcircles1p,cyl2p, bcircles2p,cyl3p,bcircles3p])
        if np.isin('STA', CART):
            fig.add_trace(STA_VMM)
            fig.add_trace(STA_LOM)
        if np.isin('VIA', CART):
            for i in viass['LINE_ID'].unique():
                via=viass[viass['LINE_ID']==i]
                fig.add_trace(go.Scatter3d(x=via['X'], y=via['Y'], z=via['Z'],
                                            hovertemplate=via['GLOBALID'],
                                            mode='lines',
                                            name='Vias y carreteras',line=dict(color='yellow',width=4),showlegend=False),)
        if np.isin('POZO', PETRO):
            fig.add_trace(Pozos)
        if len(INY)>1:
            for i in iny['CAMPO']:
                    inyc=iny[iny['CAMPO']==i]
                    fig.add_trace(go.Scatter3d(x=[float(inyc['X'])]*2, y=[float(inyc['Y'])]*2, z=[0,-1*float(inyc['prof'])],
                                hovertemplate=inyc['CAMPO'].apply(lambda x:str(x))+'<br>'
                                        'Pozos:'+inyc['POZOS'].apply(lambda x:str(x))+'<br>'
                                        'BBL:'+inyc['TOTAL_bbl'].apply(lambda x:str(x)),mode='lines',name='Inyección BBL',line=dict(color=inyc[INY],width=20,colorscale='Jet',cmax=((iny[INY])).max(),

                                cmin=((iny[INY])).min()),showlegend=False),)
            fig.update_layout(
                scene=dict(
                annotations=inyec))
        if np.isin('REZ', PETRO):
            fig.add_trace(rez)
        if np.isin('POB', CART):
                fig.add_trace(Poblaciones)
                fig.update_layout(
                scene=dict(
                annotations=Pobl))
        if np.isin('FIELD', PETRO):
            for i in campet['LINE_ID'].unique():
                f=campet[campet['LINE_ID']==i]
                attr=campet_1[campet_1['LINE_ID']==i]
                nom='Compañia:'+np.array(attr['Compañia'])[0]+'<br>Estado:'+np.array(attr['Estado'])[0]+'<br>Información:'+str(np.array(attr['INFO'])[0])
                try:
                    tip='Campo petrolífero '+np.array(attr['Campo'])[0]
                except:
                    tip='_'
                fig.add_trace(go.Scatter3d(x=f['X'], y=f['Y'], z=f['Z'],
                                hovertemplate=nom,
                                mode='lines',
                                name=tip,line=dict(color='black',width=3),showlegend=False),)
        if np.isin('LIN', PETRO):
            for i in linsis['LINE_ID'].unique():
                f=linsis[linsis['LINE_ID']==i]
                attr=linsis_1[linsis_1['LINE_ID']==i]
                try:
                    nom='ssGmStNm:'+np.array(attr['ssGmStNm'])[0]+'<br>Proyecto:'+np.array(attr['project'])[0]
                except:
                    nom='_'
                try:
                    tip=np.array(attr['owtype'])[0]
                except:
                    tip='_'
                fig.add_trace(go.Scatter3d(x=f['X'], y=f['Y'], z=f['Z'],
                                hovertemplate=nom,
                                mode='lines',
                                name=tip,line=dict(color='blue',width=3),showlegend=False),)
        if np.isin('FALL', GEOL):
            for i in faults['LINE_ID'].unique():
                fall=faults[faults['LINE_ID']==i]
                fig.add_trace(go.Scatter3d(x=fall['X'], y=fall['Y'], z=fall['Z'],
                                            hovertemplate=fall['NombreFall'],
                                            mode='lines',
                                            name='Fallas geologicas',line=dict(color='red',width=4),showlegend=False),)
        if np.isin('REAL', GEOL):
                fig.add_trace(Real)
        if np.isin('COL', GEOL):
                fig.add_trace(Colorado)
        if np.isin('MUG', GEOL):
                fig.add_trace(Mugrosa)
        if np.isin('CHO', GEOL):
                fig.add_trace(Chorros)
        if np.isin('EOC', GEOL):
                fig.add_trace(Eoceno)
        # if np.isin('PER', CART):
        #         fig.add_trace(profile_plane(x0,y0,x1,y1))
        if np.isin('SEIS', PETRO):
                fig.add_trace(SISMICA)
        fig.update_layout(autosize=False,
                        width=850, height=850,
                        margin=dict(l=50, r=50, b=50, t=50),)
        fig.update_layout(
        scene = dict(aspectratio=dict(x=2.69,y=6.24,z=(42000/693264)*EXG),
                xaxis = dict(title='Longitud(°)',nticks=10, range=[loi,los]),
                yaxis = dict(title='Latitud(°)',nticks=10, range=[lai,las],),
                zaxis = dict(title='Elevación(msnm)',nticks=10, range=[-32000,10000],),),)
        fig.update_traces(showlegend=False)
        buffer = io.StringIO()
        fig.write_html(buffer)
        html_bytes = buffer.getvalue().encode()
        encoded = b64encode(html_bytes).decode()
        return "data:text/html;base64," + encoded,START_DATE


if __name__ == "__main__":
    app.run_server(debug=True)