from urllib import *
import urllib
import numpy as np
#filters out contaminants, reverse and once identified proteins
def remove_cont_and_nan(proteins):
    proteins_filtered = proteins[(proteins['Potential contaminant'] != '+') &
                        (proteins['Only identified by site'] != '+') &
                        (proteins['Reverse'] != '+') &
                        (~(np.isnan(proteins['Ratio H/L'])))]
    return proteins_filtered

url = 'https://www.uniprot.org/uploadlists/'
def ids_to_gns(query_list, fr = 'ACC', to = 'GENENAME'):
    result_d = {}
    params = {
    'from':fr,
    'to':to,
    'format':'tab',
    'query': ' '.join(query_list)
    }
    
    data = urllib.urlencode(params)
    #request = urllib2.Request(url, data)
    request = urllib.Request(url, data)
    contact = "" #e-mail
    request.add_header('User-Agent', 'Python %s' % contact)
    #response = urllib2.urlopen(request)
    response = urllib.urlopen(request)
    page = response.read(200000)
    result = page.split('\n')
    for i in result[1:-1]:
        elements = i.split('\t')
        uid = elements[0]
        gn = elements[1]
        if uid not in result_d.keys():
            result_d[uid] = [gn]
        else:
            result_d[uid].append(gn)
    return result_d

def get_UI(iden):
    url = 'https://www.uniprot.org/uniprot/'+iden+'.txt'
    response = urllib2.urlopen(url)
    html = response.read()
    gn_start = html.find('GN   Name=')+10
    org_start = html.find('OS   ')+5
    organism = html[org_start:org_start+(html[org_start:].find('.\n'))]
    gene_name = html[gn_start:gn_start+(html[gn_start:].find(';\n'))]
    return gene_name, organism

def adjust_spines(ax, spines):
    '''adjusts spines example from matplotlib'''
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


