#para abrir los archivos
import io
#para ordenar el diccionario por valor
import operator
#para obtener log en base 2
import numpy as np
import math


#Consigna
#El programa debe descomprimir y comprimir un archivo de texto.
#Comprimir con el metodo de Fano y Markov de Orden 1
#MAYUSCULAS, puntos, spaces, newLines

def burrows_wheelerT(bwt):
    longitud = len(bwt)
    #armamos una lista donde cada item esta formado por el movimiento circular de los caracteres
    #es decir cada caracter movido a la derecha una posicion por cada iteracion
    #luego esta es ordenada alfabeticamente
    lista = sorted([bwt[i:longitud]+bwt[0:i]for i in range(longitud)])
    #posicion de la cadena original en la lista
    posicion = lista.index(bwt)
    #ahora tomamos el ultimo caracter de cada elemento de la lista para formar la cadena final
    cadena = ''
    for i in lista:
        #con -1 indicamos la ultima posicion de cada cadena
        cadena = cadena + i[-1]
    #retornamos la posicion donde se encuentra la cadena original y la cadena formada por la ultima columna luego de ordenar alfabeticamente
    return posicion,cadena

def markov_1(cadena_bwt):
    #armamos una matriz de transición
    #diccionario con las probabilidades de las transiciones
    diccionario = {}
    longitud = len(cadena_bwt)
    for i in range(longitud-1):
        if(not cadena_bwt[i] in diccionario):
            #lo que tenemos aqui es un diccionar por cada letra que a su vez tiene un diccionario con las transiciones a las siguientes letras
            #1ra vez que aparece la letra
            diccionario[cadena_bwt[i]] = {}
            d = {cadena_bwt[i+1] : 1}
            diccionario[cadena_bwt[i]].update(d)
        elif (not diccionario[cadena_bwt[i]].get(cadena_bwt[i+1])):
            #2da vez que aparece la letra pero no el siguiente caracter
            d = {cadena_bwt[i+1] : 1}
            diccionario[cadena_bwt[i]].update(d)
        else:
            diccionario[cadena_bwt[i]][cadena_bwt[i+1]]+=1
        #if(not cadena_bwt[i]+cadena_bwt[i+1] in diccionario):
        #    diccionario[cadena_bwt[i]+cadena_bwt[i+1]]=1
        #else:
        #    diccionario[cadena_bwt[i]+cadena_bwt[i+1]]+=1
    #ahora tenemos en cuenta la ultima transicion desde el ultimo simbolo hacia el inicial
    if(not cadena_bwt[-1] in diccionario):
        diccionario[cadena_bwt[-1]] = {}
        d = {cadena_bwt[0] : 1}
        diccionario[cadena_bwt[-1]].update(d)
    elif (not diccionario[cadena_bwt[-1]].get(cadena_bwt[0])):
        d = {cadena_bwt[0] : 1}
        diccionario[cadena_bwt[-1]].update(d)
    else:
        diccionario[cadena_bwt[-1]][cadena_bwt[0]]+=1
    #print(f"Diccionario: {diccionario}")
    return diccionario

def fano_mas_simbolos(transiciones, lista, clave):
    if(len(lista)==2):
        transiciones[clave][lista[0][0]] = transiciones[clave][lista[0][0]]+'0'       
        transiciones[clave][lista[1][0]] = transiciones[clave][lista[1][0]]+'1'
    elif(len(lista)==1):
        if(transiciones[clave][lista[0][0]]==''):
            transiciones[clave][lista[0][0]] = transiciones[clave][lista[0][0]]+'0'
    else:
        total=0
        for i,j in lista: 
            #calculamos el total de los caracteres actuales para hallar el 50%
            total += j
        acumulador = 0
        #para recursion
        lista1 = []
        lista2 = []
        for i,j in lista:
            if(((acumulador+j)/total) <= 0.5):
                acumulador = acumulador + j
                transiciones[clave][i] = transiciones[clave][i]+'1'       
                lista1.append((i,j))
            else:
                transiciones[clave][i] = transiciones[clave][i]+'0'
                lista2.append((i,j))
        #hacemos recursion con la partede la lista que fue asignada con 0's y la otra parte asignada con 1's
        fano_mas_simbolos(transiciones,lista1,clave)
        fano_mas_simbolos(transiciones,lista2,clave)

def arma_header(codigos,transiciones,simbolos):
    #6bits(cant. simbolos total) simbolo 4bits(pseudoascii) cant.transiciones simbolo.transicion pseudoascii long.codigo codigo...
    #6bits para indicar la cantidad de simbolos
    header = f'{len(transiciones):0{6}b}'
    for simbolo in simbolos:
        header = header + codigos[simbolo] + f'{len(transiciones[simbolo]):0{6}b}'
        for j in transiciones[simbolo]:
            header = header + codigos[j] + f'{len(transiciones[simbolo][j]):0{3}b}' + transiciones[simbolo][j]
    
    #header = f'{len(header):0{10}b}' + header
    #print(f"Header: {header}")
    return header


def fano(cadena_bwt, transiciones,pseudo_ascii):#asginar 0's y 1's a cada mejor 50%

    longitud = len(cadena_bwt)
    #obtenemos la prob. independiente de cada simbolo
    prob_indep = {}
    for i in range(longitud):
        if(cadena_bwt[i] in prob_indep):
            prob_indep[cadena_bwt[i]] += 1
        else:
            prob_indep[cadena_bwt[i]] = 1
    #..................entropia
   
    entropia = 0
    for clave in prob_indep:
        prob = (prob_indep[clave]) / longitud
        aux=0
        for i in transiciones[clave]:
            #utilizamos el valor que sera reemplazado para hallar la entropia
            prob_transicion = transiciones[clave][i]/sum(transiciones[clave][i]for i in transiciones[clave])
            #print(prob_transicion)
            aux = aux + (prob_transicion*np.log2(prob_transicion**-1))
        entropia = entropia + (prob * aux)
    print(f"Entropia del Texto Original: {entropia}")



    #..................entropia


    #proceso donde reemplazamos las probabilidades de cada clave por su codigo    
    for clave in transiciones:
        #ordenamos las probabilidades de mayor a menor en una lista    
        lista = sorted(transiciones[clave].items(), key=operator.itemgetter(1))
        for i,_ in lista:
            #eliminamos la cuenta de esa transicion y la reemplazaremos por su codigo
            transiciones[clave][i] = ''
        fano_mas_simbolos(transiciones,lista,clave)
    #print(f"Transiciones: {transiciones}")
    #-------------------fano en si es hasta la linea superior----------------------#

    #ahora que ya tenemos la coodificacion de cada transicion, utilizamos el acceso directo al diccionario para codificar
    #Colocamos el codigo correspondiente al primer estado
    cadena_codificada = str(pseudo_ascii[cadena_bwt[0]])
    for i in range(0,longitud-1):
        cadena_codificada = cadena_codificada + transiciones[cadena_bwt[i]][cadena_bwt[i+1]]
    cadena_codificada = cadena_codificada + transiciones[cadena_bwt[-1]][cadena_bwt[0]]

    header = arma_header(pseudo_ascii,transiciones,prob_indep)

    return cadena_codificada,header



def transf_ascii(linea):
    #ahora reemplazamos las secuencias de 8bits por su codificacion en ascii
    j = 0
    i = 0
    response = ""
    #tomamos de a 5 para armar en base al diccionario ascii los caracteres
    cant_bytes = int((len(linea)/5))
    while( i < cant_bytes):
        response = response + dicc2[linea[j:j+5]]        
        j = j+5
        i = i+1
    auxiliar=""
    aux_int = 0
    #comprobamos cuantos bits 0 debemos insertar
    while (j < len(linea)):
        auxiliar = auxiliar + linea[j]
        j = j+1
        aux_int = aux_int + 1 
    #insertamos los ceros hasta completar con 5 digitos binarios (coloados a la derecha)
    while(aux_int < 5):
        auxiliar += "0"
        aux_int += 1
    #concatenamos el ultimo "byte" que fue rellenado
    response = response + dicc2[auxiliar]
    return response



print("\n")
print("\n")
print("------------------------------Bienvenido al Compresor: El Picateclas------------------------------")
#en nuestro compilador no coincide la letra Ñ con la del diccionario por lo que da error al usar la letra Ñ

dicc2={'00000':'A','00001':'B','00010':'C','00011':'D','00100':'E','00101':'F','00110':'G','00111':'H','01000':'I','01001':'J','01010':'K',
'01011':'L','01100':'M','01101':'N','01110':'Ñ','01111':'O','10000':'P','10001':'Q','10010':'R','10011':'S','10100':'T','10101':'U','10110':'V',
'10111':'W','11000':'X','11001':'Y','11010':'Z','11011':',','11100':'.','11101':' ','11110':';','11111':'\n'}
pseudo_ascii={'A':'00000','B':'00001','C':'00010','D':'00011','E':'00100','F':'00101','G':'00110','H':'00111','I':'01000','J':'01001','K':'01010',
    'L':'01011','M':'01100','N':'01101','Ñ':'01110','O':'01111','P':'10000','Q':'10001','R':'10010','S':'10011','T':'10100','U':'10101',
    'V':'10110','W':'10111','X':'11000','Y':'11001','Z':'11010',',':'11011','.':'11100',' ':'11101',';':'11110','\n':'11111'}
#abrimos el archivo  a comprimir
with io.open("a_comprimir (1).txt","r") as archivo:
    cadena=archivo.read()
    #transformacion BurrowsWheeler
    indice,cadena_bwt = burrows_wheelerT(cadena)
    print(f"Indice Bwt: {indice}")
    transiciones = markov_1(cadena_bwt)
    #print(f"Transiciones: {transiciones}")
    codificado,header = fano(cadena_bwt,transiciones,pseudo_ascii)
    #print(f"Transiciones: {transiciones}")
    #CREAR ARCHIVO: comprimido.txt
    archivo = io.open("comprimido.txt","w",encoding="utf-8")
    #estos bits sobrantes se introducen cuando necesitamos completar con ceros para formar 5bits y poder hallar el correspondiente caracter
    bits_sobrantes = 5-((13 +len(header) + len(codificado)) % 5)
    linea = f'{indice:0{10}b}'+f'{bits_sobrantes:0{3}b}'+header+codificado
    linea2 = transf_ascii(linea)
    archivo.write(linea2)
    print("Ya puede ver el archivo comprimido en: \"comprimido.txt\"")
    print("\n")