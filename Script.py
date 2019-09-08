#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:24:30 2019

@author: miguel
"""


import funciones as F
print("¿Desea ver un ejemplo del funcionamiento del análisis? Y/N")
example = input().upper()

item = ""
option = ""

if example == "Y":
    print("¿Qué caso desea ver?")
    print("1 -- Cálculo de los conjuntos mínimos de corte (cutsets)")
    print("2 -- Cálculo de los Elementary Modes (modos elementales)")
    choice = input()
    if choice == str(1):
        item = "ccnf"
        option = "Cutsets"        
    elif choice == str(2):
        item = "73779"
        option = "EM"
elif example == "N":
    print("Introduzca el nombre de la proteína o el código de la reacción que desea analizar")
    print("*La proteína ha de estar contenida en la base de datos UniProt")
    print("El código ha de ser el referente a una ruta metabólica de la base de datos de Reactome")
    item = input().upper()
    if str.isdigit(item):
        item = str(item)
        option = "EM"
    else:
        option = "Cutsets"

z = F.graphs(item,option)
#F.create_html_file(z,"A",0)
