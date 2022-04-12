import pandas as pd
import plotly.graph_objects as go
import numpy as np
import skimage.io as sio
import math
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#Cilindros extracted from https://community.plotly.com/t/basic-3d-cylinders/27990/3
def cylinder(r, h,x0,y0,a =0, nt=600, nv=50):
    """
    r=radio
    h=altura
    a=altura sobre eje x y y (z)
    nt=Calidad en un eje
    nv=Calidad e otro eje
    x=coordenada x
    y=coordenada y
    """
    theta = np.linspace(0, 2*np.pi, nt)
    v = np.linspace(a, a+h, nv )
    theta, v = np.meshgrid(theta, v)
    x = r*np.cos(theta)+x0
    y = r*np.sin(theta)+y0
    z = v
    return x, y, z

def boundary_circle(r, h,x0,y0, nt=600):
    """
    r - radio del borde
    h - altura sobre el eje xy (z)
    nt-calidad de un eje
    x0-coordenada x
    y0-coordenada y
    """
    theta = np.linspace(0, 2*np.pi, nt)
    x= r*np.cos(theta)+x0
    y = r*np.sin(theta)+y0
    z = h*np.ones(theta.shape)
    return x, y, z

#Perfil geologico
def geologic_profile(x0,y0,x1,y1,url,name,colr):
    df_geo=pd.read_csv(url,delimiter=';',usecols=[2,3,4],decimal=',')
    df_geo.columns = ['X', 'Y','Z']
    x=[x0,x1]
    y=[y0,y1]
    if x0!=x1:
        slope, intercept = np.polyfit(x,y,1)
        dist=np.sqrt(((x1-x0)**2)+(y1-y0)**2)
        x_lin=np.linspace(x0,x1,int(dist*100))
        y_lin=(slope*x_lin)+intercept
    else:
        slope, intercept = np.polyfit(y,x,1)
        dist=np.sqrt(((x1-x0)**2)+(y1-y0)**2)
        y_lin=np.linspace(y0,y1,int(dist*100))
        x_lin=(slope*y_lin)+intercept
    z=[]
    d=[]
    for xln,yln in zip(x_lin,y_lin):
            df_topo_1=df_geo[(df_geo['X']<xln+0.01)&(df_geo['X']>xln-0.01)&(df_geo['Y']>yln-0.01)&(df_geo['Y']<yln+0.01)]
            dist_2=1
            forms=len(df_topo_1)
            if forms>0:
                for xt,yt,zt in zip(df_topo_1['X'],df_topo_1['Y'],df_topo_1['Z']):
                    d_prov=np.sqrt((xt-xln)**2+(yt-yln)**2)
                    if d_prov<dist_2:
                        dist_1=np.sqrt((xt-x0)**2+(yt-y0)**2)
                        dist_2=d_prov
                        z_p=zt
                if dist_2==1:
                    dist_1=np.nan
                    z_p=np.nan
            else:
                dist_1=np.nan
                z_p=np.nan
            z.append(z_p)
            d.append(dist_1)
    fig_1 = go.Scatter(x=d, y=z,
                    mode='lines',
                    name=name,
                    line=dict(color=colr, width=2))
    return fig_1
#Imagen sismica
def img_3d(url,x0,y0,x1,y1,z0,z1):
    image = sio.imread (url)
    zs,xys,_=image.shape
    zs,xys=int(zs),int(xys)
    yy = np.linspace(y0,y1, xys)
    zz = np.linspace(z1,z0, zs)
    yy, zz = np.meshgrid(yy, zz)
    x_data=list(np.linspace(x0,x1, xys))
    x_data=[x_data]*zs
    xx=np.concatenate(x_data).reshape(yy.shape)
    img = image[:,:, 1]
    ima_surface=go.Surface(x=xx, y=yy, z=zz, surfacecolor= np.flipud(img), colorscale='greys', showscale=False)
    return ima_surface
#Plano de perfil de corte
def profile_plane(x0,y0,x1,y1):
    yy = np.linspace(y0,y1, 3)
    zz = np.linspace(-15000,4000, 3)
    yy, zz = np.meshgrid(yy, zz)
    x_data=list(np.linspace(x0,x1, 3))
    x_data=[x_data]*3
    xx=np.concatenate(x_data).reshape(yy.shape)
    ima_surface=go.Surface(x=xx, y=yy, z=zz, colorscale=['red','red'], showscale=False,opacity=0.5,name='Perfil')
    return ima_surface


def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy
#Funcion para perfiles sismicos
def profile(x1,x2,y1,y2,df_sismos_1):
    x=[x1,x2]
    y=[y1,y2]
    if x1==x2:
        slope, intercept = np.polyfit(y,x,1)
        dist=np.sqrt(((x2-x1)**2)+(y2-y1)**2)
    else:
        slope, intercept = np.polyfit(x,y,1)
        dist=np.sqrt(((x2-x1)**2)+(y2-y1)**2)
    dx,dy=(x1-x2,y1-y2)
    try:
        rot=np.arctan(dx/dy)
    except:
        rot=0
    o1,o2=((x1,y1),(x2,y2))
    p1_s=rotate(o1, (x1+0.1,y1), -rot)
    p1_i=rotate(o1, (x1-0.1,y1), -rot)
    p2_s=rotate(o2, (x2+0.1,y2), -rot)
    p2_i=rotate(o2, (x2-0.1,y2), -rot)
    polygon = Polygon([p1_s, p1_i, p2_i, p2_s])
    id_ls=[]
    for x,y,id in zip(df_sismos_1[ 'LONGITUD (°)'],df_sismos_1[ 'LATITUD (°)'],df_sismos_1['Unnamed: 0']):
        if polygon.contains(Point(x,y)):
            id_ls.append(id)
    df_profile=df_sismos_1[np.isin(df_sismos_1['Unnamed: 0'], id_ls)]
    df_profile['DIST']=np.sqrt(((x1-df_profile['LONGITUD (°)'])**2)+((y1-df_profile['LATITUD (°)'])**2)) 
    return df_profile,dist
#Funcion para perfiles topograficos
def topo_profile(x0,x1,y0,y1,df_topo):
    x=[x0,x1]
    y=[y0,y1]
    if x0!=x1:
        slope, intercept = np.polyfit(x,y,1)
        dist=np.sqrt(((x1-x0)**2)+(y1-y0)**2)
        x_lin=np.linspace(x0,x1,int(dist*100))
        y_lin=(slope*x_lin)+intercept
    else:
        slope, intercept = np.polyfit(y,x,1)
        dist=np.sqrt(((x1-x0)**2)+(y1-y0)**2)
        y_lin=np.linspace(y0,y1,int(dist*100))
        x_lin=(slope*y_lin)+intercept
    z=[]
    d=[]
    for xln,yln in zip(x_lin,y_lin):
        df_topo_1=df_topo[(df_topo[0]<xln+0.01)&(df_topo[0]>xln-0.01)&(df_topo[1]>yln-0.01)&(df_topo[1]<yln+0.01)]
        dist_2=10
        for xt,yt,zt in zip(df_topo_1[0],df_topo_1[1],df_topo_1[2]):
            d_prov=np.sqrt((xt-xln)**2+(yt-yln)**2)
            # print(d_prov)
            if d_prov<dist_2:
                dist_1=np.sqrt((xt-x0)**2+(yt-y0)**2)
                dist_2=d_prov
                z_p=zt
        z.append(z_p)
        d.append(dist_1)
    fig_1 = go.Scatter(x=d, y=z,
                mode='lines',
                name='Topografía',
                line=dict(color='black', width=2))
    return fig_1

#Definicion para manejar poligonos y lineas
def lin_list(url):
    df=pd.read_csv(url)
    ls_x=[];ls_y=[];ls_z=[]
    df['Z']=df['Z']+15 #Para que sobresalgan un poco de la topografia
    for i in np.arange(0,df['LINE_ID'].max()+1):
        df_t=df[df['LINE_ID']==i]
        ls_x.append(df_t['X'].tolist())
        ls_y.append(df_t['Y'].tolist())
        ls_z.append(df_t['Z'].tolist())
    return ls_x,ls_y,ls_z

def geology(url,color_min,color_max,name_geo):
    df_geo=pd.read_csv(url,delimiter=';',usecols=[1,2,3],decimal=',')
    df_geo.columns = ['Z', 'X', 'Y']
    df_geo=df_geo[df_geo['Z']<0]
    mesh_geo=df_geo.pivot(index='Y', columns='X',values='Z')
    geology=go.Surface(z=mesh_geo.values,showscale=False, x=mesh_geo.columns, y=mesh_geo.index,showlegend=False,opacity=0.9,colorscale=[color_max,color_min],name=name_geo)
    return geology

#Funcion para graficar geologia superficial 
def geology_super(url,color,name_geo,text,TOPO):
    df_geo=pd.read_csv(url)
    df_geo=df_geo.drop_duplicates(subset=['X','Y'])
    mesh_geo=df_geo.pivot(index='Y', columns='X',values='Z')
    geology=go.Surface(z=mesh_geo.values+5,showscale=False, x=mesh_geo.columns, y=mesh_geo.index,showlegend=False,opacity=TOPO,colorscale=[color,color],name=name_geo,
                                hovertemplate=text,lighting=dict(ambient=0.2))
    return geology

def orientation(x0,y0,x1,y1):
    if x0==x1:
        if y0<y1:
            p1='S';p2='N'
        else:
            p1='N';p2='S'
    elif y0==y1:
        if x0<x1:
            p1='W';p2='E'
        else:
            p1='E';p2='W'
    else:
        alpha=np.arctan((np.abs(y0-y1))/(np.abs(x0-x1)))
        alpha=np.degrees(alpha)
        if x0<x1 and y0<y1:
            if alpha==45:
                p1='SW';p2='NE'
            elif alpha<45:
                p1='WSW';p2='ENE'
            else:
                p1='SSW';p2='NNE'
        if x0>x1 and y0>y1:
            if alpha==45:
                p1='NE';p2='SW'
            elif alpha<45:
                p1='ENE';p2='WSW'
            else:
                p1='NNE';p2='SSW'           
        if x0>x1 and y0<y1:
            if alpha==45:
                p1='SE';p2='NW'
            elif alpha<45:
                p1='ESE';p2='WNW'
            else:
                p1='SSE';p2='NNW'    
        if x0<x1 and y0>y1:
            if alpha==45:
                p1='NW';p2='SE'
            elif alpha<45:
                p1='WNW';p2='ESE'
            else:
                p1='NNW';p2='SSE'  
            
    return p1,p2

def text_scatter(SEISMO,df_sismos_1):
        ls_txt=[]
        try:
            if np.isin('LOC', SEISMO):
                ls_txt.append('Longitud:'+df_sismos_1['LONGITUD (°)'].apply(lambda x:str(x))+'°'+'<br>'
                                'Latitud:'+df_sismos_1['LATITUD (°)'].apply(lambda x:str(x))+'°'+'<br>'+
                                'Profundidad :'+df_sismos_1['PROF. (m)'].apply(lambda x:str(x))+'m <br>')
            if np.isin('MAG', SEISMO):
                ls_txt.append('Fecha :'+df_sismos_1['FECHA - HORA UTC'].apply(lambda x:str(x))+'<br>')
            if np.isin('FEC', SEISMO):
                ls_txt.append('Magnitud:'+df_sismos_1['MAGNITUD'].apply(lambda x:str(x))+'<br>'+
                'Tipo de magnitud:'+df_sismos_1['TIPO MAGNITUD']+'<br>')
            if np.isin('RMS', SEISMO):
                ls_txt.append('RMS (Segundos):'+df_sismos_1['RMS (Seg)'].apply(lambda x:str(x))+'s <br>')
            if np.isin('ERR', SEISMO):
                ls_txt.append('Error en la latitud (m):'+df_sismos_1['ERROR LATITUD (Km)'].apply(lambda x:str(x*1000))+'m <br>'+
                'Error en la longitud (m):'+df_sismos_1['ERROR LONGITUD (Km)'].apply(lambda x:str(x*1000))+'m <br>'+
                'Error en la profundidad (m):'+df_sismos_1['ERROR PROFUNDIDAD (Km)'].apply(lambda x:str(x*1000))+'m <br>')
            if len(ls_txt)==0:
                text=' '
            else:
                text=''
                for i in ls_txt:
                    text=text+i
        except:
            text='No hay sismos'
        return text