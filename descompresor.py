import io
import numpy as np

#Transformacion Burrows-Wheeler Inversa
def tbw_inversa(cadena, indice):
    longitud = len(cadena)
    #usada como matriz
    lista = [cadena[i]for i in range(longitud)]
    #la lista es un tipo de dato mutable
    lista = sorted(lista)
    #-1 por la asignacion anterior
    for _ in range(longitud-1):
        for j in range(longitud):
            lista[j] = cadena[j]+lista[j]
        lista = sorted(lista)
    print(f"Texto Original: {lista[indice]}")
    return(lista[indice])

def entropia(cadena):
    #Aquí obtenemos las transiciones del archivo comprimido
    diccionario = {}
    longitud = len(cadena)
    for i in range(longitud-1):
        if(not cadena[i] in diccionario):
            #lo que tenemos aqui es un diccionar por cada letra que a su vez tiene un diccionario con las transiciones a las siguientes letras
            #1ra vez que aparece la letra
            diccionario[cadena[i]] = {}
            d = {cadena[i+1] : 1}
            diccionario[cadena[i]].update(d)
        elif (not diccionario[cadena[i]].get(cadena[i+1])):
            #2da vez que aparece la letra pero no el siguiente caracter
            d = {cadena[i+1] : 1}
            diccionario[cadena[i]].update(d)
        else:
            diccionario[cadena[i]][cadena[i+1]]+=1
    #ahora tenemos en cuenta la ultima transicion desde el ultimo simbolo hacia el inicial
    if(not cadena[-1] in diccionario):
        diccionario[cadena[-1]] = {}
        d = {cadena[0] : 1}
        diccionario[cadena[-1]].update(d)
    elif (not diccionario[cadena[-1]].get(cadena[0])):
        d = {cadena[0] : 1}
        diccionario[cadena[-1]].update(d)
    else:
        diccionario[cadena[-1]][cadena[0]]+=1
    
    #Aqui hallamos las probabilidades independientes de cada simbolo (solo su aparicion) 
    prob_indep = {}
    for i in range(longitud):
        if(not cadena[i] in prob_indep):
            #lo que tenemos aqui es un diccionario por cada letra que a su vez tiene un diccionario con las transiciones a las siguientes letras
            #1ra vez que aparece la letra
            prob_indep[cadena[i]] = 1
        else:
            prob_indep[cadena[i]] +=1
    
    #calculamos ahora si la entropia en base a la prob de cada simbolo y sus transiciones
    entropia = 0
    for clave in prob_indep:
        prob = (prob_indep[clave]) / longitud
        aux=0
        for i in diccionario[clave]:
            #utilizamos el valor que sera reemplazado para hallar la entropia
            prob_transicion = diccionario[clave][i]/sum(diccionario[clave][i]for i in diccionario[clave])
            #print(prob_transicion)
            aux = aux + (prob_transicion*np.log2(prob_transicion**-1))
        entropia = entropia + (prob * aux)
    print(f"Entropia del Texto Comprimido: {entropia}")


def ascii_a_bits(comprimido):
    concat_bits = ""
    entropia(comprimido)
    for i in range(len(comprimido)):
        concat_bits = concat_bits + pseudo_ascii[comprimido[i]]
    return concat_bits

def datos_header(linea):
    indice = int(bytes.decode(linea[0:10].encode('utf-8')),2)
    bits_sobrantes = int(bytes.decode(linea[10:13].encode('utf-8')),2)
    cant_simbolos = int(bytes.decode(linea[13:19].encode('utf-8')),2)
    print(f'Indice Bwt: {indice}')
    transiciones = {}
    bit_actual = 19
    for _ in range(cant_simbolos):
        #NO HACERLO ENTERO SINO UNA SUCESION DE 0's y 1's (STRING) para el diccionario
        codigo_ascii = linea[bit_actual:bit_actual+5]
        #comenzamos a armar el diccionario con las transiciones entre los estados
        bit_actual += 5
        transiciones[dicc2[codigo_ascii]] = {}
        cant_transiciones = int(bytes.decode(linea[bit_actual:bit_actual+6].encode('utf-8')),2)
        bit_actual += 6
        #print("Cant. Transiciones: {0}".format(cant_transiciones))
        #print("Transiciones:")
        for _ in range(cant_transiciones):
            #String sucesion de bits para el diccionario
            ascii_siguiente = linea[bit_actual:bit_actual+5]
            bit_actual += 5
            long_codigo = int(bytes.decode(linea[bit_actual:bit_actual+3].encode('utf-8')),2)
            bit_actual += 3
            codigo_transicion = linea[bit_actual:bit_actual+long_codigo]
            bit_actual += long_codigo
            d={codigo_transicion:dicc2[ascii_siguiente]}
            transiciones[dicc2[codigo_ascii]].update(d)
            #print(transiciones)
            #print(f"Longitud de codigo: {long_codigo}")
    #leer de a bits e ir comparando con las transiciones
    datos = linea[bit_actual:-1]+linea[-1]
    descomprimida = dicc2[datos[0:5]]
    estado_actual = dicc2[datos[0:5]]
    long_datos = len(datos)#-bits_sobrantes
    i=6
    bits = datos[5]
    while(i<=long_datos-(1+bits_sobrantes)):
        #while(i<=long_datos-(bits_sobrantes+1)):  si usamos el de ascii
        if(bits in transiciones[estado_actual]):
            estado_actual=transiciones[estado_actual][bits]
            descomprimida = descomprimida + estado_actual
            bits = datos[i]
        else:
            bits = bits + datos[i]
        i += 1
    #print(f"Cadena descomprimida: {descomprimida}")
    #print(f'Transiciones: {transiciones}')
    original = tbw_inversa(descomprimida,indice)
    return original

dicc2={'00000':'A','00001':'B','00010':'C','00011':'D','00100':'E','00101':'F','00110':'G','00111':'H','01000':'I','01001':'J','01010':'K',
'01011':'L','01100':'M','01101':'N','01110':'Ñ','01111':'O','10000':'P','10001':'Q','10010':'R','10011':'S','10100':'T','10101':'U','10110':'V',
'10111':'W','11000':'X','11001':'Y','11010':'Z','11011':',','11100':'.','11101':' ','11110':';','11111':'\n'}
pseudo_ascii={'A':'00000','B':'00001','C':'00010','D':'00011','E':'00100','F':'00101','G':'00110','H':'00111','I':'01000','J':'01001','K':'01010',
    'L':'01011','M':'01100','N':'01101','Ñ':'01110','O':'01111','P':'10000','Q':'10001','R':'10010','S':'10011','T':'10100','U':'10101',
    'V':'10110','W':'10111','X':'11000','Y':'11001','Z':'11010',',':'11011','.':'11100',' ':'11101',';':'11110','\n':'11111'}
print("------------------------------Bienvenido al DEScompresor: El Picateclas------------------------------")
# 10 bits indice bwt
# 3 bits bits_sobrantes del utlimo ascii
# 6 bits cantidad de simbolos
# 5 bits pseudo ascii
# 6 cantidad de transiciones
# 5 bits pseudo ascii del codigo siguiente
# 3 bits longitud de codigo
# maximo 6 bits del codigo de transicion (indicado por los 3 bits anteriores)
with io.open("comprimido.txt","r",encoding="utf-8") as archivo:
    linea = archivo.read()
    cadena_bits = ascii_a_bits(linea)
    descomprimido = datos_header(cadena_bits)
    archivo_original = open("descomprimido11.txt","w")
    archivo_original.write(descomprimido)
    print("Ya puede ver el archivo DEScomprimido en: \"descomprimido11.txt\"")
    print("\n")