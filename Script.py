#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:24:30 2019

@author: miguel
"""
import os
import shutil
import ast
import webbrowser
import funciones as F
import networkx as nx
import matplotlib.pyplot as plt
from bioservices import Reactome

#Protein = "ccnf"

#option = "Cutsets"
#print(x)

pathway = "983168"
option = "EM"
#print(x[0],len(x[0]),"\n",x[1],len(x[1]))

def graphs(item, option):
    directory = os.path.isdir("./Graphs")
    if directory:
        shutil.rmtree('./Graphs', ignore_errors=True)
        os.mkdir("./Graphs")
    else:
        os.mkdir("./Graphs")
    
    pathways = []
    if option == "EM":
        pathways = [item]
        x = F.Elementary_Modes(item, cut = False)
    
    elif option == "Cutsets":
        x = F.cutsets(item)
        for j in x:
            pathways.append(j[2])
    
    
    for pathway in pathways:
        R = Reactome()
        
        data = R.data_discover(pathway)
        title = data['name']+" - "+pathway
        description = data['description']
        
        subdirectory = "./Data/"+pathway
    
        with open(subdirectory+"/participants.txt") as f:
            metabolites_list=ast.literal_eval(f.read())
        f.close()
        
    #        with open(subdirectory+"/reactions.txt") as f:
    #            reactions_list=ast.literal_eval(f.read())
    #        f.close()
        
        G = nx.DiGraph()
        G.clear()
        participants_index = {}
        for part in metabolites_list:
            participants_index[part[1]] = part[0]
            G.add_node(part[0])
    
    
        if option == "EM":
            EM_counter = 1
    #        print(x[0])
            for EM in x[0]:
                elemental_part_and_react = x[0][EM]
    #            print(elemental_part_and_react)
                color_nodes = []
                
                legend_nodes = {}
                for node in G:
                    colored = False
                    for P_I in participants_index:
                        if node == participants_index[P_I]:
                            if P_I in elemental_part_and_react:
                                colored = True
                            
                                color_nodes.append('cyan')
                                legend_nodes[node] = node+": "+P_I
                    if colored == False:
                        color_nodes.append('greenyellow')
                        
                legend_edges = {}
                counter = 1
                reaction_index = {}
                for relation in x[1]:
                    for i in range(len(relation)):
                        if i < len(relation)-2:
                            if relation[i] in participants_index.keys():
                                if relation[i+1] not in reaction_index:
                                    reaction_index[relation[i+1]] = str(counter)
                                    counter = counter + 1
                                G.add_edge(participants_index[relation[i]],
                                           participants_index[relation[i+2]],
                                           state=reaction_index[relation[i+1]])
                                if str(relation[i]) in str(elemental_part_and_react) and str(relation[i+2]) in str(elemental_part_and_react):
    #                                print("YES")
                                    G[participants_index[relation[i]]][participants_index[relation[i+2]]]['color'] = 'red'
                                    legend_edges["R"+str(reaction_index[relation[i+1]])] = "R"+str(reaction_index[relation[i+1]])+": "+relation[i+1]
                                else:
                                    G[participants_index[relation[i]]][participants_index[relation[i+2]]]['color'] = 'lightgray'
              
    
                LG = len(G)           
                if LG < 40:
                    node_size = 400
                    font_size = 9
                else:
                    node_size = 150
                    font_size = 6
                
                LN = len(legend_nodes)
                if LN < 11:
                    ncol_nodes = 1
                else:
                    ncol_nodes = int(str(LN)[0])+1
    
                LE = len(legend_edges)
                if LE < 11:
                    ncol_edges = 1
                else:
                    ncol_edges = int(str(LE)[0])+1
    
    #            print(len(participants_index))
    #            print(len(color_map))
                edge_color_list = [G[e[0]][e[1]]['color'] for e in G.edges() ]
                edge_labels = nx.get_edge_attributes(G,'state')
    #                fig = plt.figure()
                pos = nx.shell_layout(G)
                fig, ax = plt.subplots()
                colors=range(len(edge_labels))
                
                nx.draw(G, pos, with_labels=True,
                        node_color = color_nodes,
                        edge_color=edge_color_list,
                        node_size=node_size,font_size=font_size)
    
                nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,
                                             font_size=7)
    #            plt.savefig(directory+"/Draws/Elementary_Mode"+str(EM_counter)+".png", dpi=1000)
    
                leg = plt.legend(legend_edges.values(),bbox_to_anchor=(1.04,-0.1),
                           loc="lower left",ncol=ncol_edges,
                           fontsize=6,title = "Elementary mode reactions:\n")
        
                for j in range(LE):
                    leg.legendHandles[j].set_color('red')
        
                ax.add_artist(leg)
                Legend = plt.legend(legend_nodes.values(), bbox_to_anchor=(1.04,1),
                           loc="upper left",ncol=ncol_nodes,
                           fontsize=6, title = "Metabolites involved:\n")
    
                for k in range(LN):
                    Legend.legendHandles[k].set_color('cyan')
                plt.title(title)
    
                plt.savefig("./Graphs/Elementary_Mode"+str(EM_counter)+".png",
                            dpi=750,bbox_inches="tight")
                plt.close()
                
                EM_counter = EM_counter + 1
    #            break
                
                
        elif option == "Cutsets":
    #        print(len(participants_index))
            file = os.path.isfile(subdirectory+"List_modes.txt")
            if not file:
                F.Elementary_Modes(pathway,cut = "True")
            
            with open(subdirectory+"/List_modes.txt") as f:
                List_modes=ast.literal_eval(f.read())
            f.close()
            
            cutsets_reactions = []
            cutsets_participants = []
            
            
            for cuts in x:
                if cuts[2] == pathway:
                    for reac in cuts[0]:
                        cutsets_reactions.append(reac)
                    for met in cuts[1]:
                        cutsets_participants.append(met)
            
            print(cutsets_participants)
            color_nodes = []
            legend_nodes = {}
            for node in G:
                for P_I in participants_index:
                    if participants_index[P_I] == node:
                        
                        if P_I in cutsets_participants:
    #                        print(node,P_I)
                            color_nodes.append('cyan')
                            legend_nodes[node] = node+": "+P_I
                        else:
                            color_nodes.append('greenyellow') 
            
            counter = 1
            reaction_index = {}
            
            legend_edges = {}
            for relation in List_modes:
                for i in range(len(relation)):
                    if i < len(relation)-2:
                        if relation[i] in participants_index.keys():
                            if relation[i+1] not in reaction_index:
                                reaction_index[relation[i+1]] = str(counter)
                                counter = counter + 1
                            if relation[i+1] in cutsets_reactions:
    #                            print("YES")
                                G.add_edge(participants_index[relation[i]],
                                           participants_index[relation[i+2]],
                                           state=reaction_index[relation[i+1]]) 
                                G[participants_index[relation[i]]][participants_index[relation[i+2]]]['color'] = 'red'
                                legend_edges["R"+str(reaction_index[relation[i+1]])] = "R"+str(reaction_index[relation[i+1]])+": "+relation[i+1]
                            else:
                                G.add_edge(participants_index[relation[i]],
                                           participants_index[relation[i+2]],
                                           state=reaction_index[relation[i+1]])                               
                                G[participants_index[relation[i]]][participants_index[relation[i+2]]]['color'] = 'lightgray'
            LG = len(G)           
            if LG < 40:
                node_size = 300
                font_size = 7
            else:
                node_size = 100
                font_size = 4
            
            LN = len(legend_nodes)
            if LN < 11:
                ncol_nodes = 1
            else:
                ncol_nodes = int(str(LN)[0])
            
            LE = len(legend_edges)
            if LE < 11:
                ncol_edges = 1
            else:
                ncol_edges = int(str(LE)[0])
            edge_color_list = [G[e[0]][e[1]]['color'] for e in G.edges() ]
            edge_labels = nx.get_edge_attributes(G,'state')
    
            pos = nx.circular_layout(G)
            fig, ax = plt.subplots()
            
            nx.draw_networkx(G, pos, with_labels=True,
                    node_color = color_nodes,
                    edge_color = edge_color_list,
                    node_size=node_size,font_size=font_size)
    
            nx.draw_networkx_edge_labels(G, pos, 
                                         edge_labels=edge_labels,
                                         font_size=7)
    #        patch_nodes = mpatches.Patch(color=color_nodes, label=G.nodes)
            print(legend_edges)
            
            leg = plt.legend(legend_edges.values(),bbox_to_anchor=(1.04,-0.1),
                       loc="lower left",ncol=ncol_edges,
                       fontsize=6,title = "Minimal Cutsets:\n")
    
            for j in range(LE):
                leg.legendHandles[j].set_color('red')
                leg.legendHandles[j].set_label("A")
    
            ax.add_artist(leg)
            Legend = plt.legend(legend_nodes.values(), bbox_to_anchor=(1.04,1),
                       loc="upper left",ncol=ncol_nodes,
                       fontsize=6, title = "Metabolites involved:\n")
            for k in range(LN):
                Legend.legendHandles[k].set_color('cyan')
    #        plt.tight_layout(rect=[0,0,0.75,1])
            plt.title(title)
            plt.savefig("./Graphs/"+pathway+".png", dpi=750,
                         bbox_inches="tight")
            plt.close()
    if option == "EM":
        return([pathways, title])
    if option == "Cutsets":
        return(x)
        
def create_html_file(pathway, title, outfile):
    path = "./Graphs"
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if '.png' in file:
                files.append(os.path.join(r, file))
    for f in files:
        print(f)
    plot_file_name = files[0]
    number_of_files = len(files)-1
    counter_file = 0
#    with open(subdirectory+"/participants.txt") as f:
#        metabolites_list=ast.literal_eval(f.read())
#    f.close()
#
#    with open(subdirectory+"/reactions.txt") as f:
#        reactions_list=ast.literal_eval(f.read())
#    f.close()    
    
    from flask import Flask, render_template
    app = Flask(__name__)
    
    @app.route('/')
    def index():
        return render_template('A.html')
    
    @app.route('/crawl')
    def crawl():
        create_html_file(z[0], z[1],"A")
        return
    
    with open(f"{outfile}.html", "w") as html_file:
        html_file.write(f"\
            <!doctype html>\
            <html>\
                <head>\
                    <meta charset=\"utf-8\">\
                    <title>{title}</title>\
                    <style>\
                        body {{background-color: #f9f9f9; font-family: Helvetica, Sans-Serif;}}\
                        a {{color: blue; text-decoration: none;}}\
                    </style>\
                </head>\
                \
                <body>\
                    <h1 style=\"text-align: center;\">{title}</h1>\
                    <a target=\"_blank\" href=\"https://reactome.org/PathwayBrowser/#/R-HSA-73779\"><p style=\"text-align:center;\">Link to Reactome</a>\
                    <div style=\"display: flex;\">\
                        <object type=\"image/svg+xml\"data=\"{plot_file_name}\" height=\"520\"></object>\
                    </div>\
                    <form>\
                        <input onclick=\"window.location.href=\"/my-link/\" type=\"button\" value=\"Previous graph\" style=\"position: absolute; margin-left:40%;\"/>\
                    </form>\
                    <form>\
                        <input onclick=\"window.location.href='{files[counter_file+1]}'\" type=\"button\" value=\"Next graph\" style=\"position: absolute; margin-left:60%;\"/>\
                    </form>\
                    <a href= \'/crawl/'\ >Click me</a>\
                </body>\
            </html>")
#        if metabolites_list:
#            print("A")
#
#                    html_file.write(f"\
#                                    <h2>Metabolites:</h2>\
#                                    <ul>")
#
#                    for part in metabolites_list:
#                        name = part[0]+": "+part[1]
#                        html_file.write(f"<li><strong>{name}</strong>\</li>")
##                            <a target=\"_blank\" href=\"{go_item['xlink']['href']}\">{go_item['label']}</a></li>")
##                    if kegg:
##                        html_file.write(f"\
##                                        </ul>")
#        if reactions_list:
#
#                    html_file.write(f"\
#                                        <h2>KEGG:</h2>\
#                                        <ul>")
#
#                    for react in reactions_list:
#                        name = react[0]+": "+react[1]
#                        html_file.write(f"<li><strong>{name}</strong>\</li>")
#                            <a target=\"_blank\" href=\"{kegg_item['xlink']['href']}\">{kegg_item['label']}</a></li>")

#        html_file.write(f"\
#                    </div>\
#                </body>\
#            </html>")

pathway = "73856"
option = "EM"
z = graphs(pathway,option)
create_html_file(z[0], z[1],"A")
webbrowser.open(f"A.html")