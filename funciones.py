#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:24:10 2019

@author: miguel
"""

import os
import re
import ast
import pandas as pd
import numpy as np
from bioservices import UniProt
from bioservices import Reactome
from bioservices import ReactomeOld

def Merge(dict1, dict2): 
    res = {**dict1, **dict2} 
    return res

def pathways_Uniprot(protein):
    U = UniProt()
    R = Reactome()
    org = "organism:9606"
    rev = "reviewed:yes"
    
    
    #Busca en el archivo que contiene los nombres clave de las proteínas,
    #y extrae de ese nombre clave que ids existen en uniprot relacionados
    #con humanos. Sólo obtiene los que han sido revisados.
    
    IDs = []
    entry = protein
    x = U.search(query = protein+"+"+org+"+"+rev, columns="id,entry name,genes")
    if x == 400:
        return("No results found in UniProt")
#    print(x)
    x = x.split("\n")
    for i in x:
        if re.search("Entry",i):
            continue
        elif len(i) == 0:
            continue
        elif re.search(entry,i):
            IDs.append(i)
    IDs = set(IDs)
    dic_IDS = {}
    
    #De la lista de IDs extraída, hace una búsqueda en reactome de forma
    #que extrae los identificadores de los pathways relacionados con dichos
    #IDs.
    
    contador = 0
    for ids in IDs:
        y = U.search(ids,columns='database(Reactome)')
        y = y.split("\n")
        z1,z2,z3 = ids.split("\t")
        for j in y:
            if re.search("Cross",j):
                continue
            elif len(j) > 1:
                if z1 not in dic_IDS:
                    contador = contador + 1
                dic_IDS[z1] = j
            
        pathways = dic_IDS[z1]
        
        lista = pathways.split(";")
        l = len(lista)
        if len(lista[l-1]) < 1:
            lista = lista[:-1]
    
    dic_temp = {}
    
    for H in lista:
        h = H[6:]
        x = R.data_pathway_containedEvents(h)
        List = []
        Length = len(x)
#        print(x)
        for k in range(0,Length):
            List.append(x[k]["dbId"])
        dic_temp[h] = List
    
    list_pathways = []
    
    for key in dic_temp.keys():
        Check = False
        
        for val in dic_temp.keys():
            k = int(key)
            if k in dic_temp[val]:
                Check = True
                continue
                
        if Check == False:
            list_pathways.append(key)
    return list_pathways

def matrix_reactome(list_reactions):
    directory = os.path.isdir('./Data')
    if not directory:
        os.mkdir('./Data')
#    print("YES")
    R = Reactome()
    r = ReactomeOld()
    
#    higher_reactions = {}
#    for i in list_reactions:
#        reaction = R.data_entity_componentOf(i)
#        higher_reactions[i] = reaction[0]["stIds"]
##    print(higher_reactions)
#    ensembled_reactions = {}
#    for H in higher_reactions:
#        x = H
#        for h in higher_reactions:
#            if higher_reactions[h] == higher_reactions[H]:
#                if h !=H :
#                    x = x + "-" + h
#        if x != H:
#            ensembled_reactions[str(higher_reactions[H])[8:-2]] = x.split("-")
#        elif x == H:
#            ensembled_reactions[x] = [x]
#    ensembled_reactions = ensembled_reactions[0].split("-")
#    print(ensembled_reactions)
    
#    for ENS in ensembled_reactions:
#        ensembled_count = 0
#        list_reactions = ensembled_reactions[ENS]
##        print(list_reactions)
#        L = len(list_reactions)
#        fusioned_complexes = {}
#        fusioned_participants = {}
#        fusioned_reactions = {}
    for li in list_reactions:
#            print(li)
        subdirectory = os.path.isdir("./Data/"+li)
        if not subdirectory:
            os.mkdir("./Data/"+li)
        if os.path.isfile("./Data/"+li+"/participants.txt"):
            continue

        complexes = r.pathway_complexes(li)
        dic_complexes = {}
        Len_Comp = len(complexes)
        
        subunits = []
        com = ""
#            print("YES",li)
        for i in range(0,Len_Comp):
            if complexes[i]["schemaClass"] == "Complex":
                dic_complexes[com] = subunits
                com = complexes[i]["displayName"].upper()
                subunits = []
            
            else:
                if complexes[i]["displayName"]:
                    subunits.append(complexes[i]["displayName"].upper())
            if i == (Len_Comp - 1):
                dic_complexes[com] = subunits
        if "" in dic_complexes.keys():
            del(dic_complexes[""])
    #    keys = dic_complexes.keys()
#    
#    print (dic_complexes["ESCRT-II [cytosol]"])
            
#            print (dic_complexes)
            
#Introduce en un diccionario la lista de los participantes principales
#de la reacción. Las claves se refieren al nombre del elemento, y
#como valor se guarda una lista que contiene el identificador y 
#y el tipo de elemento.
            
        part = r.pathway_participants(li)
        
        Len_part = len(part)
        dic_participants = {}
        
        
        for j in range(0, Len_part):
            if part[j]["dbId"] == 76577:
#                    print(part[j])
                subname1,subname2 = part[j]["displayName"].split("[")
                dic_participants["AMP ["+subname2.upper()] = [part[j]["dbId"],part[j]["schemaClass"].upper()]
            elif part[j]["dbId"] == 164121:
#                    print(part[j])
                subname1,subname2 = part[j]["displayName"].split("[")
                dic_participants["AMP ["+subname2.upper()] = [part[j]["dbId"],part[j]["schemaClass"].upper()]
            elif part[j]["displayName"] == "Ub [nucleoplasm]":
                dic_participants[part[j]["displayName"].upper()] = ["68524",part[j]["schemaClass"].upper()]
            else:
                dic_participants[part[j]["displayName"].upper()] = [part[j]["dbId"],part[j]["schemaClass"].upper()]
        
#        print ("\n\n",dic_participants,len(dic_participants))
    
    #Se guarda en un diccionario las reacciones que participan de la reacción
    #en cuestión, cuya clave es el id y el valor el nombre.
        
        react = R.data_pathway_containedEvents(li)
        
#        subpaths = {}
        dic_react = {}
        Len_react = len(react)    
        
        for k in range(0,Len_react):
    #        print(R.data_pathway_containedEvents(react[k]["dbId"]))
#            if R.data_pathway_containedEvents(react[k]["dbId"]) == 404:
            dic_react[react[k]["dbId"]] = react[k]["displayName"].upper()
#            else:
#                subpaths[react[k]["dbId"]] = react[k]["displayName"]
            
#        print("\n\n",dic_react)
#        print("\n\n\n\n\n",dic_reactions)
    
#            print("\n\n",len(dic_react))
        
#        Continue = False    
##            print(L)
#        if L == 1:
#            Continue = True
#            fusioned_complexes = dic_complexes
#            fusioned_participants = dic_participants
#            fusioned_reactions = dic_react
#        else:
#            ensembled_count = ensembled_count + 1
#            if ensembled_count == 1:
#                fusioned_complexes = dic_complexes
#                fusioned_participants = dic_participants
#                fusioned_reactions = dic_react
#            else:
#                fusioned_complexes = Merge(dic_complexes,fusioned_complexes)
#                fusioned_participants = Merge(dic_participants,fusioned_participants)
#                fusioned_reactions = Merge(dic_react,fusioned_reactions)
#                if ensembled_count == L:
#                    Continue = True
#                    li = ENS
#                    subdirectory = os.path.isdir("./Data/"+li)
#                    if not subdirectory:
#                        os.mkdir("./Data/"+li)           
                    
    
    #A partir de cada uno de los participantes se extraen las reacciones
    #directas sobre las que participa. Lo usual es que exista una reacción 
    #de la que dicho participante es un reactivo (INPUT) y otra de la que es
    #un producto (OUTPUT). Esto se guarda en un diccionario denominado
    #"dic_relations"
#        if Continue == True:
        dic_relations = {}
            
        for l in dic_participants.keys():
            x = dic_participants[l][0]
            reaction = R.data_entity_componentOf(x)
            
            if reaction == 404:
                continue
    
            OP = ""
            IP = ""
            for m in reaction:
                if m["type"] == "output":
                    OP = m["names"]
                elif m["type"] == "input":
                    IP = m["names"]
    
            IP_list = []
            OP_list = []
    
            if IP:
                IP = list(set(IP))
                for ip in IP:
                    ip = ip.upper()
                    if ip in dic_react.values():
                        IP_list.append(ip)
            if OP:
                OP = list(set(OP))
                for op in OP:
                    op = op.upper()
                    if op in dic_react.values():
                        OP_list.append(op)
            if (IP_list and OP_list):
                dic_relations[l] = [{"SUBSTRATUM OF":IP_list},{"PRODUCT OF":OP_list}]
                
            elif (IP_list and not OP_list):
                dic_relations[l] = {"SUBSTRATUM OF":IP_list}
            
            else:
                dic_relations[l] = {"PRODUCT OF":OP_list}
        
#        print("\n\n",dic_relations)
        
    #Se generan dos listas que contienen los nombres de las reacciones y de los
    #participantes
    
        reactions = dic_react.values()
        reactions = list(reactions)
        
        participants = dic_participants.keys()
        participants = list(participants)
    
    #    print("\n\n",reactions[0],"\n",len(reactions))
    #    print("\n\n",participants,"\n",len(participants))
    
    #Para cada reacción se extrae el código html que presenta la página web de
    #Reactome. Esto se hace con la finalidad de extraer los elementos catalíticos,
    #y lo que regulan positiva y negativamente a cada reacción. Se guardan para
    #ello trozos de texto de código html, para posteriormente analizarlo.
    #Además, se analiza la presencia de participantes de la reacción que pudieran
    #realizar la función de sustrato o producto en más de una ocasión. Es decir,
    #que para dicha reacción se necesiten 2 o más moléculas de un mismo tipo de
    #participante (metabolito).
    
        count = 0
        ids_reactions = []
        dic_elements = {}
        check = []
        list_repeated = []
        
        for RE in dic_react.keys():
            if os.path.isdir('./Temp') == False:
                os.mkdir('./Temp')
            
            index = "R"+str(count)
            ids_reactions.append([index,dic_react[RE]])
            count = count + 1
    #        print(count)
    #        print("\n\n",RE)
    
            html = r.bioservices_get_reactants_from_reaction_identifier(RE)
            
            f = open("./Temp/html.txt","w+")
            f.write(str(html))
            f.close() 
            
            Cat = ""
            Positive = ""
            Negative = ""
            Requirements = ""
            
            values = []
            sections = []
            
            for line in open("./Temp/html.txt"):
    #            print(line)
                contador = 0
                for i in line:
#                    print(i,"\n\n")
                    if line[contador:contador+8] == "Catalyst":
                        Cat = line[contador:contador+2000]
                    elif line[contador:contador+10] == "Positively":
                        Positive = line[contador:contador+2000]
                    elif line[contador:contador+10] == "Negatively":
                        Negative = line[contador:contador+2000]
                    elif line[contador:contador+12] == "Requirements":
                        Requirements = line[contador:contador+2000]                     
                    elif line[contador:contador+3] == " x ":
                        values.append(line[contador-1])
                        sections.append(line[contador:contador+500])
                    contador = contador + 1
    #        print(RE)
#            print(Requirements)
                
            temp_catalytic = []
            temp_positive = []
            temp_negative = []
            temp_requirements = []
            
            if Cat:
                t = len(Cat)
            elif Positive:
                t = len(Positive)
            else:
                t = len(Negative)
            
            for P in participants:
                edited_P = False
                for letter in P:
                    if letter == "'":
                        edited_P = True
                        L1,L2 = P.split("'")
                if edited_P == True:
#                    print(L1)
#                    print(L2)
                    temp_P = L1+"\\'"+L2
#                    print(temp_P)
                    pat = len(temp_P)
#                    print(len(temp_P))
                else:
                    pat = len(P)
                    temp_P = P
                for k in range(0,len(sections)):
                    T = len(sections[k])
                    for pos in range(0, T - (pat - 1)):
                        if temp_P == sections[k][pos:pos+pat].upper():
                            list_repeated.append([values[k],P,dic_react[RE]])
    
                
    #            if P == "Glc [nucleoplasm]":
    #                if RE == 170810:
    #                    print("YES")
                for position in range(0, t - (pat - 1)):
    #                print(pos)
            #        print(Cat[pos:pos+p])
    
                    if temp_P == Cat[position:position+pat].upper():
                        temp_catalytic.append(P)
                        check.append(P)
    
                    elif temp_P == Positive[position:position+pat].upper():
                        temp_positive.append(P)
                        check.append(P)
    
                    elif temp_P == Negative[position:position+pat].upper():
                        temp_negative.append(P)
                        check.append(P)
                    if Requirements:
                        if temp_P == Requirements[position:position+pat].upper():
                            temp_requirements.append(P)
   
#            print(temp_requirements)
            if temp_requirements:
                for req in temp_requirements:
                    temp_dic_relations = []
                    if req in dic_relations:
#                        print(dic_relations[req])
                        if len(dic_relations[req]) == 1:
                            if "SUBSTRATUM OF" in dic_relations[req]:
                                temp_dic_relations = dic_relations[req]["SUBSTRATUM OF"]
                                del(dic_relations[req])
                                temp_dic_relations.append(dic_react[RE])
                                dic_relations[req] = {"SUBSTRATUM OF":temp_dic_relations}
                            else:
                                temp_dic_relations = {}
                                temp_dic_relations = dic_relations[req]["PRODUCT OF"]
                                del(dic_relations[req])
                                dic_relations[req] = [{"SUBSTRATUM OF":dic_react[RE]},{"PRODUCT OF":temp_dic_relations}]
                        else:
                            temp_dic_relations_S = dic_relations[req][0]["SUBSTRATUM OF"]
                            temp_dic_relations_P = dic_relations[req][0]["PRODUCT OF"]
                            del(dic_relations[req])
                            temp_dic_relations_S.append(dic_react[RE])
                            dic_relations[req] = [{"SUBSTRATUM OF":temp_dic_relations_S},{"PRODUCT OF":temp_dic_relations_P}]                    
                    else:
                        dic_relations[req] = {"SUBSTRATUM OF":dic_react[RE]}
    
            if (temp_catalytic or temp_positive or temp_negative):
                dic_elements[dic_react[RE]] = [{"Catalytic":list(set(temp_catalytic))},
                            {"Positive":list(set(temp_positive))},
                            {"Negative":list(set(temp_negative))}]
    #    print(dic_elements)
        
    #Para cada participante se extrae la información sobre las reacciones en
    #las que participa para construir la matriz estequiométrica. Las columnas de
    #la matriz representan las reacciones, mientras que las filas lo hacen con
    #los metabolitos
    #
    #Si el participante es consumido en una reacción, en la matriz se escribirá 
    #un -1. Por el contrario, si es producido se escribirá un +1. En caso de no 
    #existir relación entre un parcipante y una reacción dados, se escribirá un 0.
    #Antes de colocar los casos donde hay consumo o producción, se comprueba que
    #el metabolito no interacciona de forma múltiple (en una reacción participan
    #varios metabolitos del mismo tipo). Si es así, se sutituye el número 1 por
    #el referenciado en "dic_repeated".
    
        check = list(set(check))
    #    print(len(check))
        
        ids_participants = []
        matrix = []
        contador = 0
        
#        print("\n",list_repeated,"\n")
#        print("\n",dic_repeated,"\n")
#        print(check)
        for p in participants:
            index = "M"+str(contador)
#            if p in dic_relations.keys():
            ids_participants.append([index,p])
            contador = contador + 1
#            print("\n\n",p)
    
            temp_input = []
            temp_output = []
            line = []
        #    print(dic_relations[p],"\n",p)
            
    #        print(dic_elements.values())
    
            if (p in check) and (p not in dic_relations.keys()):
                for q in reactions:
                    line.append(0)
                matrix.append(line)
                continue
            elif len(dic_relations[p]) == 1:
                if "SUBSTRATUM OF" in dic_relations[p]:
                    temp_input = dic_relations[p]["SUBSTRATUM OF"]
                elif "PRODUCT OF" in dic_relations[p]:
                    temp_output = dic_relations[p]["PRODUCT OF"]
            
            elif len(dic_relations[p]) > 1:
                temp_input = dic_relations[p][0]["SUBSTRATUM OF"]
                temp_output = dic_relations[p][1]["PRODUCT OF"]
        
        
            for q in reactions:
                I = False
                O = False        
    
                if temp_input:
        #            print(temp_input)
                    for t_i in temp_input:
                        if t_i == q:
                            I = True
        
                if temp_output:
        #            print(temp_output)
                    for t_o in temp_output:
                        if t_o == q:
                            O = True            
                
                Repeated = False
                for rep in list_repeated:
                    if (rep[2] == q and rep[1] == p):
                        if I == True:
                            number = -int(rep[0])
                        elif O == True:
                            number = +int(rep[0])
                        elif (I == True and O == True):
                            number = 0
                        else:
                            number = 0
                        line.append(number)
                        Repeated = True
    
                if Repeated == True:
                    continue
                else:
                    if I == True and O == True:
                        line.append(0)
                    elif I == True and O == False:
                        line.append(-1)
                    elif I == False and O == True:
                        line.append(+1)
                    elif I == False and O == False:
                        line.append(0)
            
        #    print("\n",p,"\n",line)
            matrix.append(line)
        
                        
        #print("\n",dic_catalytic)
        f = open("./Data/"+li+"/complexes.txt","w+")
        f.write(str(dic_complexes))
        f.close() 
        
        f = open("./Data/"+li+"/participants.txt","w+")
        f.write(str(ids_participants))
        f.close() 
        
        f = open("./Data/"+li+"/reactions.txt","w+")
        f.write(str(ids_reactions)+"\n")
        f.close()
        
        f = open("./Data/"+li+"/relations.txt","w+")
        f.write(str(dic_relations)+"\n")
        f.close()        
        
        path = "./Data/"+li+"/matriz.txt"
        df = pd.DataFrame(matrix)
        df.to_csv (path,
                   header = False,
                   index = False,
                   sep = ' ')
    
        columns = "\t"
        f = open("./Data/"+li+"/matriz2.txt","w+")
        for REACT in ids_reactions:
            columns = columns + " " + REACT[0]
        f.write(columns+"\n")
        count_rows = 0
        for Line in matrix:
            rows = ""
            space = ""
            count_columns = 0
            for item in Line:
                if count_columns > 9:
                    space = "   "
                else:
                    space = "  "
                if len(str(item)) == 1:
                    rows = rows + space + str(item)
                else:
                    rows = rows + space[:-1] + str(item)
                count_columns = count_columns + 1
            rows = ids_participants[count_rows][0] + "\t" + rows + "\n"
            f.write(rows)
            count_rows = count_rows + 1
        f.close()
    return(list_reactions)
            #    print(ids_participants,"\n",ids_reactions)
            
            
def cutsets(word):
    word = word.upper()
    reactions_Uniprot = pathways_Uniprot(word)
    matrix_reactome(reactions_Uniprot)
    cutsets = []
    for pathway in reactions_Uniprot:
#        pathway = "674695"
        directory = "./Data/"+pathway

        with open(directory+"/complexes.txt") as f:
            complexes_dic=ast.literal_eval(f.read())
        f.close()

        with open(directory+"/participants.txt") as f:
            metabolites_list=ast.literal_eval(f.read())
        f.close()

        with open(directory+"/reactions.txt") as f:
            reactions_list=ast.literal_eval(f.read())
        f.close()

        with open(directory+"/relations.txt") as f:
            relations_dic=ast.literal_eval(f.read())
        f.close()
         
        complexes_involved = []
        
        L = len(word)
        for comp in complexes_dic:
            for subunit in complexes_dic[comp]:
                if subunit[0:L] == word:
                    complexes_involved.append(comp)
        
        participants_involved = []
        for part in metabolites_list:
            if part[1][0:L] == word or part[1] in complexes_involved:
                participants_involved.append(part[1])

        reactions_related = {}

        A = np.loadtxt("./Data/"+pathway+"/matriz.txt")
        
        #    Este bucle genera un diccionario que introduce como clave
        #    una reacción, y como valor la o las reacciones de las que
        #    depende
        for participant in participants_involved:
            for i in range(A.shape[0]):
                part = metabolites_list[i][1]
                if part == participant:
                    temp_react = []
                    elemental_react = []
                    for j in range(A.shape[1]):
                        react = reactions_list[j][1]
                        if int(A[i][j]) < 0:
                            temp_react.append(react)
                        elif int(A[i][j]) > 0:
                            elemental_react.append(react)
                    if temp_react and elemental_react:
                        for elem in elemental_react:
                            if elem in reactions_related.keys():
                                add = reactions_related[elem]+temp_react
                                add = list(set(add))
                                del(reactions_related[elem])
                                reactions_related[str(elem)] = add
                            else:
                                reactions_related[str(elem)] = temp_react
                    if temp_react and not elemental_react:
                        reactions_related[part] = ["S",temp_react]
                    elif elemental_react and not temp_react:
                        reactions_related[part] = ["P",elemental_react]
#        print(participants_involved)
#        print(reactions_related,"\n")
            
        for key in reactions_related:
#            print("YES")
            introduce = True
            for KEY in reactions_related:
                if reactions_related[KEY][0] == "S":
                    if reactions_related[KEY][1] == key:
                        introduce = False
                if key in reactions_related[KEY]:
                    introduce = False
            if introduce == True:
                for P in relations_dic:
                    if P in participants_involved:
                        if len(relations_dic[P]) == 1:
                            
                            if "PRODUCT OF" in relations_dic[P]:
                                
                                for react in relations_dic[P]["PRODUCT OF"]:
#                                    print([P,react])
                                    if P == key:
                                        
                                        cutsets.append([react,P,pathway])
                        else:
                            if "PRODUCT OF" in relations_dic[P][1]:
                                for react in relations_dic[P][1]["PRODUCT OF"]:
                                    if react == key:
                                        cutsets.append([react,P,pathway])
#                    cutsets.append([key,pathway])
#                    print(key,"\n\n")
        
        if len(reactions_related) == 1:
            if len(participants_involved) == 1:
                if participants_involved[0] in reactions_related.keys():
                    Participant = participants_involved[0]
                    if reactions_related[Participant][0] == "P":
                        for R in reactions_related[Participant][1]:
                            cutsets.append(R,participant,pathway)
                    elif reactions_related[Participant][0] == "S":
                        for r in reactions_related[Participant][1]:
                            cutsets.append([r,participant,pathway])

#    print(cutsets,"\n\n")
    temp_cutsets = cutsets
    cutsets = []
    for CUT in temp_cutsets:
        metabolites = []
        reactions = []
        pathway = []
        for cut in temp_cutsets:
            if cut[-1] == CUT[-1]:
                if cut[0] not in metabolites: 
                    metabolites.append(cut[0])
                if cut[1] not in reactions: 
                    reactions.append(cut[1])
        pathway = [metabolites,reactions,CUT[-1]]
        if pathway not in cutsets:
            cutsets.append(pathway)
#    print(cutsets)    
    return(cutsets)

def Elementary_Modes(Pathway,cut):
    Pathway = Pathway.upper()
    matrix_reactome([Pathway])
    directory = "./Data/"+Pathway
    
    if os.path.isfile(directory+"/Elementary_modes.txt"):
        with open(directory+"/Elementary_modes.txt") as f:
            Elementary_modes=ast.literal_eval(f.read())
        f.close()
        
        with open(directory+"/List_modes.txt") as f:
            List_modes=ast.literal_eval(f.read())
        f.close()
        
        return(Elementary_modes,List_modes)
    else:    
        with open(directory+"/participants.txt") as f:
            metabolites_list=ast.literal_eval(f.read())
        f.close()
        
        with open(directory+"/reactions.txt") as f:
            reactions_list=ast.literal_eval(f.read())
        f.close()
        
        with open(directory+"/relations.txt") as f:
            relations_dic=ast.literal_eval(f.read())
        f.close()
    
        #l = len(word)
        #for com in complexes_dic:
        #    if word == com[0:1]:
        #        for element in complexes_dic[com]:
        #            INPUT.append(element)
        #        INPUT.append(com)
        #if not INPUT:
        #    INPUT = [word]
        
        #print(INPUT)
        participants_involved = []
    
        for part in metabolites_list:
            participants_involved.append(part[1])
        
    #    print(participants_involved,"\n\n")
        
        reactions_related = {}
    
        A = np.loadtxt("./Data/"+Pathway+"/matriz.txt")
        
        #    Este bucle genera un diccionario que introduce como clave
        #    una reacción, y como valor la o las reacciones de las que
        #    depende
    #    0 = Filas
        for participant in participants_involved:
            for i in range(A.shape[0]):
                part = metabolites_list[i][1]
                if part == participant:
                    temp_react = []
                    elemental_react = []
                    for j in range(A.shape[1]):
                        react = reactions_list[j][1]
                        if int(A[i][j]) < 0:
                            temp_react.append(react)
                        elif int(A[i][j]) > 0:
                            elemental_react.append(react)
                    if temp_react and elemental_react:
                        for elem in elemental_react:
                            if elem in reactions_related.keys():
                                add = reactions_related[elem]+temp_react
                                add = list(set(add))
                                del(reactions_related[elem])
                                reactions_related[str(elem)] = add
                            else:
                                reactions_related[str(elem)] = temp_react
                    if temp_react and not elemental_react:
                        reactions_related[part] = ["S",temp_react]
                    elif elemental_react and not temp_react:
                        reactions_related[part] = ["P",elemental_react]
        
        reactions_involved = []
        for x in reactions_list:
            reactions_involved.append(x[1])
        
        circular_pathways = []
        for metabolite in relations_dic:
            if len(relations_dic[metabolite]) == 2:
        #        print("YES")
        #        print(relations_dic[metabolite][0])
        #        for subs in relations_dic[metabolite][0]["SUBSTRATUM OF"]:
        #            for METABOLITE in relations_dic:
        #                
        #                if len(relations_dic[METABOLITE]) == 2:
        #                    for SUBS in relations_dic[METABOLITE][1]["PRODUCT OF"]:
        ##                        print(SUBS)
        #                        if subs == SUBS:
        #                            if metabolite not in circular_metabolites:
        #                                circular_metabolites.append([METABOLITE,SUBS])
                for pro in relations_dic[metabolite][1]["PRODUCT OF"]:
                    for METABOLITE in relations_dic:
                        if len(relations_dic[METABOLITE]) == 2:
                            for PRO in relations_dic[METABOLITE][0]["SUBSTRATUM OF"]:
                                if pro == PRO:
                                    if metabolite not in circular_pathways:
                                        circular_pathways.append([METABOLITE,PRO])      
        #print(circular_pathways)
        
        start_modes = {}
        for start in reactions_related:
            if reactions_related[start][0] == "S":
                temp = []
                seq = reactions_related[start][1]
                temp = [seq[i:i+1] for i in range(0, len(seq), 1)]
                start_modes[start] = temp
        
#        print(start_modes)
        
        circular_metabolites = []
        circular_reactions = []
        for circular in circular_pathways:
            temp_circular = []
            t_c = circular[0]
            circular_metabolites.append(t_c)
            circular_reactions.append(circular[1])
            for CIRCULAR in circular_pathways:
                T_C = CIRCULAR[0]
                if t_c == T_C:
                    
                    if CIRCULAR[1] not in temp_circular:
        #                print(CIRCULAR[1])
                        temp_circular.append([CIRCULAR[1]])
            start_modes[t_c] = temp_circular
        
#        print(start_modes)
        #print(len(start_modes))
        
        all_modes_and_participants = {}
        some_ways = []
        temp_ways = []
        for key in start_modes:
            all_ways = []
            
            for n in start_modes[key]:
                temp_ways = [n]
                some_ways = temp_ways
                for i in range(len(participants_involved)+len(reactions_list)):
                    
                    if not some_ways:
                        continue
                    temp_ways = []
                    for j in range(len(some_ways)):
                        temp = []  
                        
                        posibilities = []
                        select = some_ways[j][-1]
        #                print(select)
            #            continue
            #            if "SUPT6H [NUCLEOPLASM]" in str(select):
            #                print("YES")
        #Si la última reacción de some_ways es igual a alguna de las
        #reacciones que produzcan uno o varios metabolitos en cuestion,
        #y esto metabolitos solo sean producidos por dicha reaccion,
        #se añade a la lista de reacciones. Esto se considera un final
        #de dicha lista.
                        for R in reactions_related:
                            if reactions_related[R][0] == "P":
                                for reaction in reactions_related[R][1]:
                                    if select == reaction:
                                        posibilities.append(R)
        #                print(posibilities)
                        BREAK = False
                        if posibilities:
                            for POS in posibilities:
                                temp = some_ways[j]
#                                if POS == key:
#                                    print(key,"\n",temp,"\n",POS)
#                                    BREAK = True
                                if POS not in temp:
                                    temp = temp + [POS]
                                if temp not in all_ways:
                                    all_ways.append(temp)
                        if BREAK == True:
                            break

                        temp = []
                        posibilities = []
                        for P in relations_dic:
                            if len(relations_dic[P]) == 1:
                                if select in reactions_involved:
                                    if "PRODUCT OF" in relations_dic[P]:
                                        for react in relations_dic[P]["PRODUCT OF"]:
                                            if react == select:
                                                posibilities.append(P)
                                elif select in participants_involved: 
                                    if "SUBSTRATUM OF" in relations_dic[P]:
                                        if P == select:
                                            for react in relations_dic[P]["SUBSTRATUM OF"]:
                                                Introduce = True
        #                                        print(some_ways[j])
                                                if key in circular_metabolites:
        #                                            print(react)
        #                                            if str(react) == "2-4 NT.BACKTRACKING OF POL II COMPLEX ON THE TEMPLATE LEADING TO ELONGATION PAUSING":
        #                                                print("YES")
                                                    if str(react) in str(circular_reactions):
        #                                                print("YES")
                                                        Introduce = False
                                                        temp = some_ways[j]
                                                        if temp not in all_ways:
                                                            all_ways.append(temp)
                                                if Introduce == True: 
                                                    posibilities.append(react)
                            else:
        #                        print("YES")
                                if select in reactions_involved:
                                    if "PRODUCT OF" in relations_dic[P][1]:
                                        for react in relations_dic[P][1]["PRODUCT OF"]:
                                            if react == select:
                                                posibilities.append(P)
                                elif select in participants_involved: 
                                    if "SUBSTRATUM OF" in relations_dic[P][0]:
                                        if P == select:
                                            for react in relations_dic[P][0]["SUBSTRATUM OF"]:
                                                Introduce = True
                                                if key in circular_metabolites:
        #                                            print(some_ways[j])
        #                                            if str(react) == "2-4 NT.BACKTRACKING OF POL II COMPLEX ON THE TEMPLATE LEADING TO ELONGATION PAUSING":
                                                                                                    
                                                    if str(react) in str(circular_reactions):
                                                    
        #                                                    print(react)
        #                                                    print(some_ways[j][-2],"\n")
                                                        Introduce = False
                                                        temp = some_ways[j]
                                                        if temp not in all_ways:
                                                            all_ways.append(temp)
                                                if Introduce == True: 
                                                    posibilities.append(react)
                        temp = []
        #                print(posibilities)
                        if posibilities:
                            for pos in posibilities:
                                temp = some_ways[j]
                                
                                if pos in temp:
                                    all_ways.append(temp)
                                    continue
                                temp = temp + [pos]
                                if temp not in temp_ways:
                                    temp_ways.append(temp)
            
                    some_ways = temp_ways
                
            all_modes_and_participants[key] = all_ways
        
        #print("\n\n\n\n\n\n\n\n\n\n",all_modes_and_participants)
        #print(len(all_modes_and_participants))
        
        List_modes = []
        for K in all_modes_and_participants:
            complete = []
            for List in all_modes_and_participants[K]:
        #            print(List)
                List = [K]+List
                complete.append(List)
            List_modes = List_modes + complete
        
        #print("\n\n\n\n\n\n\n\n\n\n")
        #for LIst in List_modes:
        #    print(LIst,"\n")
        #print("\n\n\n\n\n\n\n\n\n\n",List_modes)
        #print(List_modes)
#        print(len(List_modes))
        
        temp_LM = List_modes
        List_modes = []
        List_counter = 0
        for List in temp_LM:
            List_counter = List_counter + 1
        #    print(List_counter)
            Introduce = True
            for LIST in temp_LM:
                if List == LIST:
                    continue
                if str(List)[1:-1] in str(LIST)[1:-1]:
                    Introduce = False
            if Introduce == True:
                if str(List)[1:-1] not in str(List_modes):
                    List_modes.append(List)
#        print(List_modes)
#        print(len(List_modes))
        
        f = open(directory+"/List_modes.txt","w+")
        f.write(str(List_modes))
        f.close()
        
        if cut == "True":
            return
        
        index = []
        for P in List_modes:
            counter1 = 0
            
            for i in range(len(P)):
                if len(P) == 3:
                    substring1 = [P[0]]+[P[1]]+[P[2]]
                else:
                    if counter1 < len(P)-2:
                        substring1 = [P[counter1]]+[P[counter1+1]]+[P[counter1+2]]
                    else:
                        counter1 = counter1 + 1
                        continue
                
                for Q in List_modes:
                    counter2 = 0
                    for j in range(len(Q)):
                        if len(Q) == 3:
                            substring2 = [Q[0]]+[Q[1]]+[Q[2]]
                        else:
                            if counter2 < len(Q)-2:
                                substring2 = [Q[counter2]]+[Q[counter2+1]]+[Q[counter2+2]]
                            
                            else:
                                counter2 = counter2 + 1
                                continue
        #                if substring1[0] in participants_involved:
                        if substring1[0] == substring2[0]:
                            if substring1[1] == substring2[1]:
                                if substring1[2] != substring2[2]:
                                    if substring1 not in index:
                                        index.append(substring1)
                                    if substring2 not in index:
                                        index.append(substring2)
                        if substring1[0] != substring2[0]:
                            if substring1[1] == substring2[1]:
                                if substring1[2] == substring2[2]:
                                    if substring1 not in index:
                                        index.append(substring1)
                                    if substring2 not in index:
                                        index.append(substring2)
        #                elif substring1[1] in participants_involved:        
        #                    if substring1[0] == substring2[0]:
        #                        if substring1[1] == substring2[1]:
        #                            if substring1[2] != substring2[2]:
        #                                if substring1 not in index:
        #                                    index.append(substring1)
        #                                if substring2 not in index:
        #                                    index.append(substring2)
        #                    if substring1[0] != substring2[0]:
        #                        if substring1[1] == substring2[1]:
        #                            if substring1[2] == substring2[2]:
        #                                if substring1 not in index:
        #                                    index.append(substring1)
        #                                if substring2 not in index:
        #                                    index.append(substring2)               
                        counter2 = counter2 + 1
                counter1 = counter1 + 1
        
#        print(index)
        temp_index = index
        index = []
        join_index = []
        for I in temp_index:
            if I[1] in participants_involved:
                if I not in index:
                    index.append(I)
                continue
            start = []
            end = []
            for IN in temp_index:
                if IN[1] in participants_involved:
                    continue
                if I[1] == IN[1]:
                    if (IN[0]) not in start:
                        start.append(IN[0])
                    if (IN[2]) not in end:
                        end.append(IN[2])
            candidate = [start,I[1],end]
            if candidate not in join_index:
                join_index.append(candidate)
#        print(join_index)
#        print(index)
        
        #repeat = 0
        #central = []
        #for iN in index:
        #    if iN[1] not in central:
        #        central.append(iN[1])
        ##        print(iN[1])
        #repeat = len(central)
        ##print(repeat)
        #
        #incompatibilities = []
        #if repeat > 1:
        #    for x in it.product(index,repeat = repeat):
        #        x = list(x)
        #        Introduce = True
        #        for i in range(len(x)):
        #            for j in range(len(x)):
        #                if i == j:
        #                    continue
        #                if x[i][1] == x[j][1]:
        #                    Introduce = False
        #        if Introduce == True:
        #            if str(sorted(x)) not in str(incompatibilities):
        #                incompatibilities.append(sorted(x))
            
        
        #else:
        #    for IND in index:
        #        incompatibilities.append([IND])
        #print(incompatibilities)
        
        Elementary_modes = {}
        Elementary_modes_counter = 0
        total_counter = 0
        for R in List_modes:
            print("Entrada")
            total_counter = total_counter + 1
        
        #    print(R,"\n")
        #    continue
        #    if total_counter == 2:
        #        break
        
        #    for incomp in incompatibilities:
        #        impossible = []
        #        for INCOMP in incompatibilities:
        #            for INC in INCOMP:
        #                if str(INC)[1:-1] not in str(incomp):
        #                    impossible.append(INC)
                
                
            temp_P_T_S = [R]
#            if str(temp_P_T_S) != "[['ATP [NUCLEOPLASM]', 'UNWINDING OF DNA FOR THE NASCENT TRANSCRIPT: SECOND TRANSITION', 'PPI [NUCLEOPLASM]']]":
#                continue   
#                print("\n\n",impossible,"\n\n")
            definitive_list = []
            pathways_to_study = []
            for i in range(len(List_modes)):
                if temp_P_T_S == []:
                    break
                for pathway in temp_P_T_S:
                    introduce = True
                    for react in pathway:
                        select_ends = []
                        select_starts = []
                        for ind in join_index:
                            reaction_index = ind[1]
                            
                            if react == reaction_index:
#                                print(react,"\n")
        #                        continue
                                introduce = False
                                pre_ends = ind[0]
                                post_starts = ind[2]
                                
                                select_ends = []
                                select_starts = []
                                if len(pre_ends) > 1:
                                    for p_e in pre_ends:
                                        select_ends.append([p_e,reaction_index])
                                else:
                                    select_ends = [[pre_ends[0],reaction_index]]
                                if len(post_starts) > 1:
                                    for p_s in post_starts:
                                        select_starts.append([reaction_index,p_s])
                                else:
                                    select_starts = [[reaction_index,post_starts[0]]]
#                                    print(select_ends)
#                                    print(select_starts,"\n\n")
                        if select_ends == [] and select_starts == []:
                            continue
                        else:
                            for S in List_modes:
        #                        Continue = False
        #                        for imp in impossible:
        #                            if str(imp)[1:-1] in str(S):
        #                                Continue = True
        #                        if Continue == True:
        #                            continue
                                for j in range(len(S)):
                                    substring_S = S[j-1:j+1]
                                    for se in select_ends:
                                        if str(se)[1:-1] == str(substring_S)[1:-1]:                               
                                            if str(S[0:j])[1:-1] not in str(pathways_to_study):
                                                pathways_to_study.append(S[0:j])
                                                definitive_list.append(S[0:j])
        
                                            if str(se[0]) not in str(definitive_list):
                                                definitive_list.append([se[0]])                                      
                                    for ss in select_starts:
                                        if str(ss)[1:-1] == str(substring_S)[1:-1]:
                                            if str(S[j:len(S)])[1:-1] not in str(pathways_to_study):
                                                pathways_to_study.append(S[j:len(S)]) 
                                                definitive_list.append(S[j:len(S)])
                                            if str(ss[0]) not in str(definitive_list):
            #                                                print("SS\n",ss[0])
                                                definitive_list.append([ss[0]])                  
        #                print("PPPPPPPPP\n",pathways_to_study,"\n",len(pathways_to_study))
        #                print("\nDDDDDDDD\n",definitive_list,"\n",len(definitive_list))
        #            break
                    if introduce == True:
                        if pathway not in definitive_list:
                            definitive_list.append(pathway)
        #                    print("YES")
        
#                print("\n111111111\n",pathways_to_study,len(pathways_to_study))
        
                if pathways_to_study:
                    for path in pathways_to_study:                      
                        if str(path) == "[]":
                            pathways_to_study.remove(path)
                        elif len(path) == 1:
                            if path not in definitive_list:
                                definitive_list.append(path)
                            pathways_to_study.remove(path)
                            
#                print("\n222222222\n",pathways_to_study,len(pathways_to_study))
#                print("\n3333333333\n",temp_P_T_S,len(temp_P_T_S))
                if temp_P_T_S == pathways_to_study:
                    temp_P_T_S = []
                    continue
                
#                print("\n3333333333\n",temp_P_T_S)
        #            print(path)
        #            print(definitive_list)
        #            elif len(path) == 2:
        #                if path[0] == "Introduce":
        #                    if path not in definitive_list:
        #                        definitive_list.append(path)
        #                    pathways_to_study.remove(path)
        
        #        print(R,"\n")
        #        print(temp_P_T_S,"\n")
        #        print(pathways_to_study,"\n\n")
        #        
        #        print(len(temp_P_T_S))
        #        print(len(pathways_to_study))
                
                
                temp_P_T_S = pathways_to_study
        #        pathways_to_study = []
        
                
                
        #        pathways_to_study = []
        
        #        print(pathways_to_study)
        #        print("\n",definitive_list,"\n\n\n")
                
                
#            print("\n\n",definitive_list,len(definitive_list))        
            temp_def_list = definitive_list
            definitive_list = []
            for defin in temp_def_list:
                for element in defin:
                    if element not in definitive_list:
                        definitive_list.append(element)
#            print("\n\n",definitive_list,len(definitive_list))
        
            INTRODUCE = True
            for Keys in Elementary_modes:
        #            print(definitive_list)
                if sorted(definitive_list) == sorted(Elementary_modes[Keys]):
                    INTRODUCE = False            
            if INTRODUCE == True:
                Elementary_modes_counter = Elementary_modes_counter + 1
                Elementary_modes[Elementary_modes_counter] = definitive_list
            print(total_counter," de ",len(List_modes))
        
        #    break
    #    print(len(Elementary_modes))    
        f = open(directory+"/Elementary_modes.txt","w+")
        f.write(str(Elementary_modes))
        f.close()
        return(Elementary_modes,List_modes)
    #print("\n\n",Elementary_modes,"\n",len(Elementary_modes))             
    
    
    #print(len(index))
    

