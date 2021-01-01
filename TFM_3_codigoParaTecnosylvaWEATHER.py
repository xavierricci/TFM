#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Thu Nov 26 20:04:37 2020

@author: xaf



# RECOLECCION DE DATOS CLIMATOLOGICOS DIARIOS 


Esta funcion se ha creado para que Tecnosliva la haga correr sobre su fuente de archivos netcdf con el historico de datos climatologicos.

- IMPUTS:
    
     1. dfFuegosConInfoTopo.csv

     Un dataframe con todos los fuegos y sus fechas
     
       - fecha 	fecha de inicio del incendio (solo el anyo)
       - long_F 	coord_Igni, coordenadas(long en WGS 84) del punto de centroide del perimetro del fuego
       - lat_F 	coord_Igni, coordenadas(lat en WGS 84) del punto de centroide del perimetro del fuego
       - datos 	ya se han borrado los puntos random sin datos topograficos. ahora solo hay Con datos y Nuevas Coordenadas (que algunos datos random que se movieron un poco del punto original y se reescribieron las coordenadas) 	
       - indice 
       
     2. dfRandomNOfuegoConDatos
     Una carpeta con listados de los puntos random SIN FUEGO x anyo, organizado en un dataframe con 5 columnas:

        - fecha 	fecha de inicio del incendio (solo el anyo)
        - long_F 	coord_Igni, coordenadas(long en WGS 84) del punto de centroide del perimetro del fuego
        - lat_F 	coord_Igni, coordenadas(lat en WGS 84) del punto de centroide del perimetro del fuego
        - datos 	ya se han borrado los puntos random sin datos topograficos. ahora solo hay Con datos y Nuevas Coordenadas (que algunos datos random que se movieron un poco del punto original y se reescribieron las coordenadas) 	
        - indice 
        
     3. nc
     Una carpeta que representa los netcdfs con informacion metereologica 
      
      
      
      - _____*** HAY Q REDIRECCIONAR nc a la direccion donde estan los datos meteorologicos !!! ***_____
     
     
     
     
     4. poligonoDerecogida 
     Poligono SHAPEFILE que se va creando a cada recogida
     
     5. poligonoMeteo
     2 Poligono SHAPEFILE que limin el espacio de recogida de los datos meteorologicos 
         

- PROCESO: 

    1. A partir del listado de fuegos (de los que se ha confirmado que tenemos datos topograficos), usando sus fechas se importa la información climatológica del día (netcdf)

    2. A partir de las coordenadas se genera una matrices de fuego y una de no fuego. ambas de 27x27x4
    
    cada pixel representa 500m por lo que la matriz representa 13.5 km

    las variables metereologicas (con sus limites) son :
    
    
      - _______________________*** HAY Q CONFIRMARLAS O CAMBIARLAS!!! ***________________________
     
        - v_LFM=[0,300]
        - v_U10=[-20,+20]
        - v_V10=[-20,+20]
        - v_MO1=[0,1]
        
        
        
        
    3. Todas las matrices de fuego de un mismo anyo se gurardan en una lista que se guarda en un archivo en formato pickel, que a su vez se gurda en una carpeta llamada * matricesFMeteo
    
    4. Simultaneamente, por cada fuego con datos metereologicos, a partir de un listado de puntos random se se genera una matriz con los datos metereologicos del mismo dia en el punto random. Todas las matrices de sitios sin fuego de un mismo anyo se gurardan en una lista que se guarda en un archivo en formato pickel, que a su vez se gurda en una carpeta llamada * matricesNoFMeteo

    5. Simultaneamente, por cada fuego anazlizado se genera un df con 9 atributos:
    A los 5 iniciales anyade:
    
        - meteo_F  Que inidca si el punto original tiene datos meteorologicos 
        - LO_NoF  La longitud del punto random 
        - LA_NoF  La latidut del punto random	
        - met_NoF Que inidca si el punto random tiene datos meteorologicos 

    estos datos acumulados por anyo en df se guradan en la carpeta dfMETEOtotal
    
    (este df serivira para posteriormente localizar y combinar las matrices con datos topograficos con las de datos meterologicos)
    
    
- OUTCOME: (YA HAY CARPETAS PREPARADAS)

    matricesFMeteo

    matricesNoFMeteo

    matricesFMeteo
    


'''


#2. CREACION DE UNA MATRIZ DE 5 X 5 DE LOS PIXELS ORIGINALES DE WHEADER QUE SON DE CASI 3KM POR LADO (UNOS 15KM) 
#Y QUE POR LO TANTO RECOGEN MAS DE LOS 13.5 KM QUE NECESITAMOS
###_____________________________________


import numpy as np
import pandas as pd
from netCDF4 import Dataset
from PIL import Image
import os

import datetime
import pickle



#This will return a 5x5 matrix 
def GetRaster(lon,lat,path,varname):
    with Dataset(path, "r") as ncin:
    
       
        
        #ii,jj=LatLon_To_index_BruteForce(path,lat,lon)
        ii,jj=LatLon_TO_index(lat,lon)
        #print(f'{ii} {jj}')
      
        ras=np.zeros([5,5])
        
        myVar = ncin.variables[varname]
        
        for _i in range(5):
            for _j in range(5):
                
                ras[_i][_j]=np.mean(myVar[:,ii-2+_i , jj-2+_j])
        
        return ras

def LatLon_TO_index(lat, lon):
    nrow=320
    ncol=256
    latmin=33.89
    latmax=42.45
    lonmin=-125.42
    lonmax=-116.59
    
    # for a given lat,lon we need to know the values ii,jj in the matrix
    # to do that we compute the distance in degrees for each cell in the matrix (iiStep)
    # and then compute ii,jj
    iiStep = (latmax - latmin) / nrow
    jjStep = (lonmax - lonmin) / ncol
    ii = (
        latmax - lat
    ) / iiStep  # this is latmax-lat because I think ii starts on the N (check it)
    jj = (lon - lonmin) / jjStep  # this is lon-lonmin because jj starts on the West
    return (nrow-ii, jj)



def LatLon_To_index_BruteForce(path,lat,lon):
    with Dataset(path, "r") as ncin:
        nrow=320
        ncol=256
        XLAT = ncin.variables["XLAT"]
        XLON = ncin.variables["XLONG"]
        closest_ii=-1
        closest_jj=-1
        
        dis=999999
        for ii in range(nrow):
            for jj in range(ncol):
                _lat=XLAT[0,ii,jj]
                _lon=XLON[0,ii,jj]
               
                distance=(lat-_lat)**2+(lon-_lon)**2
                
                if distance<dis:
                  
                    dis=distance
                    closest_ii=ii
                    closest_jj=jj
    return (closest_ii,closest_jj)




def cordAmatrizRapido(longx,latx,netcdfDELDIA):
  
    atributoS =['lfm',  'U10','V10', 'mean_wtd_moisture_1hr'] 

    matriz3D=[]
    for x in atributoS:
        m=GetRaster(longx,latx,netcdfDELDIA,x) 
        matriz3D.append(m)
    return matriz3D



# puesto que estamos con fuegos de california y sabemos que el nc solo corta al este y al sur nos aseguramos de que los fuegos mas alla no enttran.
# este limite se hace mas de 6km antes del borde del netcdf para que pueda recoger datos alrededor del punto requerido
def estaDentro(long,lat): #=boundary):
    estaDentro=True
    if lat < 34 or long > -117:
        estaDentro=False
    return estaDentro




#FOLER_TO_NC='nc/'
FOLER_TO_NC="D:/0_FIRECAST/weather/pgehistorical/"
#FOLER_TO_NC="/home/computeadmin/projects/weather/pgehistorical/"

# 1. RECOGEMOS LOS DATOS CLIMATOLOGICOS DEL DIA EN NC Y LOS PASAMOS A DF


# 4 CREAMOS UNA FUNCION QUE A PARTIR DE UNA MATRIZ DADA LA NORMALIZA 
#Y LA CONVIERTE EN UNA MATRIZ DE 27X27 
#DE LA QUE POSTERIORMENTE PODREMOS EXTRAER matrices 9 matrices de 25x25 PARA LA RED NEURONAL 

def prepararMatrix(CapaMatriz, ValorMinimoPosible=-999,ValorMaximoPosible=-999, numeroPixels=27, normalizada=True):

    
   
    #pasamos el array a 1 dimension para poder cambiar Todos los datos extremos por el valor maximo
    a1D = CapaMatriz.reshape((len(CapaMatriz)*len(CapaMatriz[0])),1)

    #creamos un contador de extremos (si hubiera muchos, posiblemente habria un problema con el supuesto valor maximo)
    contadorExtremos0 = 0
    if ValorMinimoPosible != -999:
        for valor in range(len(a1D)):
            if a1D[valor] < ValorMinimoPosible:
                a1D[valor] = ValorMinimoPosible

                contadorExtremos0 = contadorExtremos0+1
       
        #printLog(f"Se han encontrado {contadorExtremos0} elementos extremos que se han substituido por el valor minimo: {ValorMinimoPosible}")       
        if contadorExtremos0>2:
            printLog(printLog(f"Se han encontrado {contadorExtremos0} elementos extremos que se han substituido por el valor minimo: {ValorMinimoPosible}"))       

    else:
        #Si estos valores no estuvieran saldrá un aviso
        ValorMinimoPosible=0
        printLog("falta el valor minimo posible de este atributo")

 
    contadorExtremos = 0

    if ValorMaximoPosible != -999:
        for valor in range(len(a1D)):
            if a1D[valor] > ValorMaximoPosible:
                a1D[valor] = ValorMaximoPosible
                contadorExtremos = contadorExtremos+1

        #printLog(f"Se han encontrado {contadorExtremos} elementos extremos que se han substituido por el valor maximo: {ValorMaximoPosible}")       
        if contadorExtremos0>2:
            printLog(printLog(f"Se han encontrado {contadorExtremos0} elementos extremos que se han substituido por el valor minimo: {ValorMinimoPosible}"))       
                  
    else:
        #Si estos valores no estuvieran saldrá un aviso
        ValorMaximoPosible=50
        printLog("falta el valor maximo posible de este atributo")
        
        

    try:

        #devolvemos el array a la matriz original SIN LOS VALORES EXTREMOS
        XArray = a1D.reshape(len(CapaMatriz),len(CapaMatriz[0]))
    #printLog(f"Sin extremos:\n{XArray}")

        #guardamos los max y minimos de la matriz original
        miniMat=XArray.min()
        maxiMat=XArray.max()
    #printLog(miniMat)

        XArrayA0 = (XArray-XArray.min()) #esto mueve el rango de 0 a max
        arrayNorm0_1 = XArrayA0/XArrayA0.max() #esto lo normaliza
    #printLog(f"Normalizada:{arrayNorm0_1}") #Esto es con lo que alimentaremos la imagen para la creacion o fusion de N pixels


        #convertimos el array en una imagen con valor min = 0 azul y max = 255 amarillo
        img = Image.fromarray((arrayNorm0_1 * 255).astype(np.uint8))
    #plt.matshow(img)
        #convertimos la imagen en matriz
        img2 = img.resize((numeroPixels, numeroPixels), Image.BICUBIC) #ANTIALIAS) #BILINEAR) NEAREST

        arrayImagen = np.array(img2, copy=True)

        #printLog(f"Nuevo Tamanyo:\n{arrayImagen}")

        # CONVERSION DE LA IMAGEN EN COLOR MODIFICADA O NO, A UNA MATRIZ EN BASE A LOS MAXIMOS Y MINIMOS DE LA MATRIZ ORIGINAL

        casiNormalizado=arrayImagen/arrayImagen.max()


        miniFoto=casiNormalizado.min()
        maxiFoto=casiNormalizado.max()

        for x in range(len(casiNormalizado)):
            for y in range(len(casiNormalizado[x])):
                if casiNormalizado[x][y] == miniFoto:
                    casiNormalizado[x][y] = miniMat
                elif casiNormalizado[x][y] == maxiFoto:
                    casiNormalizado[x][y] = maxiMat
                else: 
                    casiNormalizado[x][y] = ((casiNormalizado[x][y])*(maxiMat-miniMat)) + miniMat
            conValoresReales=casiNormalizado


        if normalizada==False :        
            return conValoresReales     

        else:
            #printLog(casiNormalizado)

            rango = ValorMaximoPosible - ValorMinimoPosible

            matrizNormalizadaEnAbsolutos = (conValoresReales-ValorMinimoPosible)/rango

            return matrizNormalizadaEnAbsolutos 
    except:
        printLog('matriz no ajustable')
        return -9999

    
    

    
    
# 5. FUNCION DE RECOGIDE DATOS 


'''
ES NECESARIO ESTABLECER/CONFIRMAR LOS MINIMOS Y MAXIMOS:
v_LFM=[0,100]
v_U10=[-20,+20]
v_V10=[-20,+20]
v_MO1=[0,1]
'''


def ncCoord2cnnMatrix3DDD(long,lat,netcdfDELDIA):
#def ncCoord2cnnMatrix3DDD(long, lat,geodataframe,boundary):
    
    if estaDentro(long, lat)==False:
        return -9999
    
    matrizLFM_27=cordAmatrizRapido(long,lat,netcdfDELDIA)
    if matrizLFM_27 == -9999:
        return -9999

 ##############________________________________________________REVISAR_________________________________########

    #### ES NECESARIO ESTABLECER LOS MINIMOS Y MAXIMOS:
    v_LFM=[0,300]
    v_U10=[-20,+20]
    v_V10=[-20,+20]    
    v_MO1=[0,1]
   
##############________________________________________________REVISAR_________________________________########


    atributoS =[v_LFM, v_U10, v_V10, v_MO1]  
    

    #ESTA MATRIZ REPRESENTA 15KM a 3KM POR PIXEL (25pixels),  PERO NOSOTROS SOLO QUEREMOS 12.5KM
    # asi que dividimos cada pixel en 6 partes: 25*6=150

    m27x27x4=[]
    for atr in range(4):
        m150=prepararMatrix(matrizLFM_27[atr],atributoS[atr][0],atributoS[atr][1],150)
  
    #       if m150==-9999:
    #           return -9999
        
#Ahora cada pixel representa 100m, cotamos 8 y 7 pixels por lado para reducir 1.5km que nos sobran
        X1 = m150[8:153, 8:153]
    #plt.matshow(X1)
    #plt.colorbar()

    #agrupamos los valores para devolver una matiz de 25 x 25 pixels que representan 500m de lado, cada pixel
    #Ahora el maximos siempre es 1 i el minimo 0 pq ya esta normalizada!!
        m27F=prepararMatrix(X1,0,1,27)
        m27x27x4.append(m27F)
    #plt.matshow(m25F)
    #plt.colorbar()
    return m27x27x4










# 6. A PARTIR DEL DF CON TODOS LOS FUEGOS CON DATOS SE PREPARA EL DF PARA LOS DATOS FINALES X ANYOS
#ESTAS SON LAS COORDENADAS DEFINITIVAS CON INFORMACION TOPOGRAFICA (ya tienen maties organizadas x anyos en este orden)

def prepararNuevoDF():
    direccion=('dfFuegosConInfoTopo.csv')
    dfFuegosTopo = pd.read_csv(direccion, index_col=0)

    dfFuegosTopo["meteo_F"]='sinDatos'

    dfFuegosTopo["indice"]=range(len(dfFuegosTopo))
    dfFuegosTopo["LO_NoF"]=0.0
    dfFuegosTopo["LA_NoF"]=0.0
    dfFuegosTopo["met_NoF"]='sinDatos'
    dfFuegosTopo=dfFuegosTopo.reset_index()
    del dfFuegosTopo["index"]
    return dfFuegosTopo
#prepararNuevoDF()



#AQUI EMPIEZA LA FUNCION DE CREACION DE MATRICES DE DATOS METEOROLOGICOS para cada anyo:
#AQUI EMPIEZA LA FUNCION DE CREACION DE MATRICES DE DATOS METEOROLOGICOS para cada anyo:
def matizMeteoXanyo(dfDelAnyo,anyo):#, boundary):
    
    matrizMeteoFuegoAnyo=[]
    matrizMeteoNO_FuegoAnyo=[]
    print(1111111111111111)  
    #abrimos la lista de puntos random sin fuego posibles del anyo a tratar
    direccion=('dfRandomNOfuegoConDatos/dfFinalNOfuego'+str(anyo)+'.csv')
    dfNOFuegoAnual = pd.read_csv(direccion, index_col=0)
    

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    #### PARA HACER PROEBAS ###

    #for posicion in range(0,5):  #PARA HACER PRUEBBAS SACAMOS SOLO LOS 5 PRIMEROS DE CADA ANYO
    #PARA LA VERSION FINAL HAY QUE BORRAR EL ANTERIOR Y ACTIVAR EL SIGUIENTE
    for posicion in range(len(dfDelAnyo)):       

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

        
        printLog(anyo)
        actual=dfDelAnyo['fecha'][posicion]
             
        printLog(f"actual: {actual}")
        
        #ncdf= FOLER_TO_NC+"wrf_LT_d03_2010-11-13.nc"
        ncdf= FOLER_TO_NC+"wrf_LT_d03_"+actual+".nc"

        if os.path.isfile(ncdf)==False:
            printLog(f'{ncdf} Does not exist') 

            continue

        
        #detectamos si la fecha se repite, si es así, no hace falta cargar los datos de nuevo
        #para evitar problemas con el primero de la lista
        if posicion==0:
            posicionB=posicion+60
            anterior=dfDelAnyo['fecha'][posicionB]
        else:
            anterior=dfDelAnyo['fecha'][posicion-1]
    
      
        
        if actual != anterior:  
            

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    

        
    # ABRIMOS LOS DATOS METEOROLOGICOS/ AQUI SE TIENEN QUE INDICAR EL CAMINO A LOS DATOS!
            #CONEXION A LOS DATOS DE TECNOSYLVA
           

#############################



#### PRUEBA CON UNA UNICA DIRECCION!!!!



            try:
                #ncdf= FOLER_TO_NC+"wrf_LT_d03_2010-11-13.nc"
                ncdf= FOLER_TO_NC+"wrf_LT_d03_"+actual+".nc"

               
                if os.path.isfile(ncdf)==False:
                    printLog(f'{ncdf} Does not exist') 
                    
                #datospoli=meteoAmatriz(datosEntrecomillas=ncdf)

        #activar la linea anterior y borrar la seiguiente
                #datospoli=meteoAmatriz(datosEntrecomillas="nc/wrf_LT_d03_2010-11-13.nc")
                printLog(f"Aqui se cargaria el ncdf de datos metereologicos diarios como datospoli: {ncdf}")    


            except Exception as e: 
                
                printLog('Failed: '+ str(e))
                printLog(f'_____No hay datos de esta fecha:{actual}')
                #datospoli= []
                continue


        else:
            #en caso de que se repita la fecha pero no haya habido informacion metereologica
            #if len(datospoli)==0:
             #   printLog(f'_____MISMA FECHA:{actual} SIN DATOS')
              #  datospoli= []
            printLog(f"Continuar usando los datos anteriores: {ncdf}")



        # SE REPITE LA FECHA Y HAY DATOS
        long=float(dfDelAnyo['long_F'][posicion])
        lat=float(dfDelAnyo['lat_F'][posicion])
        #printLog(long,lat)

         
    #aqui recogemos la matriz para esta posicion    
    
        
        matrizMeteoFuegodia=ncCoord2cnnMatrix3DDD(long,lat,ncdf)
        if matrizMeteoFuegodia==-9999:
            printLog("no tiene datos meteorologicos")

    #SI EL FUEGO NO TIENE DATOS SE ACABA AQUÍ y PASA AL SIGUIENTE FUEGO 

        else:
            #imprimimos el primer numero central de la primera matriz para ver si funciona (que no se repita el mismo)
            printLog(matrizMeteoFuegodia[0][0][14])
            #guardamos la matriz con las matrices del anyo
            matrizMeteoFuegoAnyo.append(matrizMeteoFuegodia)

            #anotamos que tiene matriz en el df
            dfDelAnyo.at[posicion,'meteo_F']="Con datos"

    #EMPIEZA LA RECOGIDA DE DATOS PARA LA MATRIZ NO FUEGO para cada fuego con datos meteo
    
            #printLog(f"Quedan {len(dfNOFuegoAnual)} puntos random del anyo {anyo} ")
                
            contadorCoordRandom=0
            
            matriz=False
            while matriz==False:   

                #recogemos las coordenadas de la primera posicion de la lista.
                longR=dfNOFuegoAnual['long_F'][contadorCoordRandom]
                latR=dfNOFuegoAnual['lat_F'][contadorCoordRandom]
                #printLog(longR,latR)
                
                
                matrizNO=ncCoord2cnnMatrix3DDD(longR,latR,ncdf)#datospoli, boundary)

                if matrizNO != -9999:

                    matriz=True
                    matrizMeteoNO_FuegoAnyo.append(matrizNO)

                    #indicamos en el df
                    dfDelAnyo.at[posicion,'met_NoF']="Con datos"
                    dfDelAnyo.at[posicion,'LO_NoF']=longR       
                    dfDelAnyo.at[posicion,'LA_NoF']=latR
                    
                    #guardar la posicion del dfNOFuegoAnual
                    dfDelAnyo['dfNOF']=contadorCoordRandom
                    
                    
                    #imprimimos el primer numero central de la primera matriz para ver si funciona (que no se repita el mismo)
                    printLog(matrizNO[0][0][14])
                    
                    #borramos las coordenadas ya usadas
                    dfNOFuegoAnual=dfNOFuegoAnual.drop(contadorCoordRandom)
                    dfNOFuegoAnual=dfNOFuegoAnual.reset_index()
                    del dfNOFuegoAnual['index']
                    
                else:
                    printLog("no hay datos meteorologicos del punto random. Pasamos al siguiente")
                    contadorCoordRandom+=1

                    
    printLog('AÑO COMPLETADO!!!')
    # GURADAMOS LAS MATRICES CON LOS DATOS ANUALES EN LAS CARPETAS
    direccion=('dfMETEOtotal/dfMETEOtotal'+str(anyo)+'.csv')
    dfDelAnyo.to_csv(direccion)

    
    direccion=('matricesFMeteo/matricesFMeteo'+str(anyo)+'.pkl')
    with open(direccion, 'wb') as f:
        pickle.dump(matrizMeteoFuegoAnyo, f)

    direccion=('matricesNoFMeteo/matricesNoFMeteo'+str(anyo)+'.pkl')
    with open(direccion, 'wb') as f:
        pickle.dump(matrizMeteoNO_FuegoAnyo, f)

    #SI SE QUIERE COMPROBAR QUE LOS DATOS SE GUARDAN CORRECTAMENTE
    #return dfDelAnyo


def printLog(a):
    print(a)
    f.write(f'{a} \n')


def main():
   
# 7. LANZAMOS LA FUNCION PARA CADA ANYO DE LA LISTA 


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    

#CAMBIAR CUANDO TENGAMOS MAS ANYOS AHORA SOLO 2019

    #anyos = list(range(2019,2020))
    anyos = list(range(2000,2020))
    
    for anyoAbierto in anyos:
        dfFuegosTopo2=prepararNuevoDF()

        
        dfFuegosTopo2["fecha"]= pd.to_datetime(dfFuegosTopo2["fecha"])

        
        date_time_obj = datetime.datetime.strptime(str(anyoAbierto), '%Y')
        date_time_obj2 = datetime.datetime.strptime(str(anyoAbierto+1), '%Y')
        dfDelAnyo=dfFuegosTopo2[(dfFuegosTopo2["fecha"]>=date_time_obj)&(dfFuegosTopo2["fecha"]<date_time_obj2)]

        
        dfDelAnyo["fecha"] = dfDelAnyo["fecha"].astype(str) 

        
        dfDelAnyo=dfDelAnyo.reset_index()
        del dfDelAnyo['index']
        dfDelAnyo['dfNOF']=0
        
        #printLog(f'Anyo: {anyoAbierto}')
        #printLog(dfDelAnyo[0:2])    

        matizMeteoXanyo(dfDelAnyo, anyoAbierto)#, boundary)

    
    
        
        
   
    
if __name__ == "__main__":
    with open("LOG.txt", "w") as f:
        main()

    

def testAnyo(anyo):  
   # 8. Comprobacion de que el sistema funciona correctamente
    direccion=('dfMETEOtotal/dfMETEOtotal'+str(anyo)+'.csv')
    meteo2002 = pd.read_csv(direccion, index_col=0)
    
    nofff=len(meteo2002[meteo2002['met_NoF']=='Con datos'])
    printLog(f'El df de 2002 tiene {nofff} puntos de NoFuego con datos')
    
    
    with open('matricesNoFMeteo/matricesNoFMeteo'+str(anyo)+'.pkl', 'rb') as f:
        ml2002 = pickle.load(f)   
    printLog(f'El mismo anyo tiene {len(ml2002)} matrices guradadas')
    
    return meteo2002 