import numpy as np
from pylab import *
from matplotlib.colors import colorConverter
from matplotlib.collections import RegularPolyCollection
from matplotlib.widgets import RadioButtons
import matplotlib.gridspec as gridspec
import math

from GeneSetRankedEnrichment23 import *

LIVER=0
accuracy=15

if LIVER:
    m = GeneSetRankedEnrichment(1, 1, 5, 500, 0, 0.9, 0.01, 3, 0.01, 
                                   0.15, 0.0001, 5, 0.12, 10, 25, 2, '3-cell-hep-d12_versus_2-cell-hep-d12_abs.txt',
                                   'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                   'set_bounds_file.txt',
                                   'results_GOEnrichment.txt', 'LIVER')
else:
    m = GeneSetRankedEnrichment(2, 2, 5, 500, 0, 0.9, 0.01, 200, 0.01, 
                                   0.15, 0.001, 2, 0.6, 1000, 1000, 2, 'bkz_gene_list.csv',
                                   'human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv',
                                   'set_bounds_file.txt',
                                   'results_GOEnrichment.txt', 'HDF')
(C,genes_in_C)=m.runMe(0,0,1,0,0,0)
unique_genes=m.loadvar('unique_genes')
lug=len(unique_genes)
T=m.loadvar('T')
#large_goids=m.loadvar('large_goids')
Tset=set(range(len(T)))
term_names=m.loadvar('term_names')
cc=[len(T[i]) for i in xrange(len(T))]
(genes_list, bounds, level_list) = m.get_sets_and_genes()
if m.s==0:
    bounds[-1]=lug
(G,_, dict_G)=m.getG(unique_genes,genes_list,T,Tset,bounds)
levels=[level_list[genes_list.index(unique_genes[i])] for i in xrange(len(unique_genes))]

plot_title=("GenGO" if m.kind==1 else "Bayesian")+r" (belief$=$"+str(m.belief)+", "
if m.kind==1:
    plot_title+=r'$p=$%s, $q=$%s, $\alpha=$%s)' % (round(m.p,4),round(m.q,4),round(m.a,3))
else:
    plot_title+=r'$\alpha=$%s, $\beta=$%s, $p=$%s)' % (round(m.alpha,4),round(m.beta,4),round(m.prob,6))
plot_title+=", $s=$%s,\n%s data, %s$\geq$size of GO terms$\geq$%s, %s terms" % (m.s,m.out_str,m.top_cutoff if m.top_cutoff>m.cutoff else r'$\infty$',m.cutoff,len(C))


if m.s==1:
    category_mat=[[0,1],[0],[1]]
elif m.s>0:
    category_mat=[range(m.s+1),range(m.s),[m.s],[0]]
else:
    category_mat=[range(len(dict_G)+1),range(len(dict_G)),[len(dict_G)]]


##################Real start of code, above all copied from GeneSetRankedEnrichment
xlabel_ind=0
ylabel_ind=3
geneslabel_ind_x=1
geneslabel_ind_y=1

keys_genes=['all genes','perturbed genes','unperturbed genes']
for i,c in zip([lug,sum([len(g) for g in G[:-1]]),len(G[-1])],range(3)):
    keys_genes[c]+=' (%s)' % i
if m.s>1:
    keys_genes.append('highest perturbed genes (%s)' % len(G[0]))

#Create Basic data
total_annotations=[cc[c] for c in C]
perturbed_mat=[[cc[c] for c in C],[len(T[c])-len(set(T[c])&G[-1]) for c in C],
               [len(set(T[c])&G[-1]) for c in C],[len(set(T[c])&G[0]) for c in C]]
mean_gene_level_mat=[m.mean_gene_level(C,G,T,levels),m.mean_gene_level(C,G,T,levels,range(m.s)),
                     m.mean_gene_level(C,G,T,levels,[m.s]),m.mean_gene_level(C,G,T,levels,[0])]
pvalues=[math.log(item,10) for item in m.p_values_of_list(C, T, G)]#np.array(res[6])
term_names=[term_names[c] for c in C]#'RNA splicing & mRNA splicing', 'respiratory electron transport chain', 'mitotic prometaphase', 'translational initiation', 'S phase', 'DNA-dependent transcription', 'regulation of Rho protein signal transduction', 'RNA localization', 'anaphase-promoting complex-dependent proteasomal ubiquitin-dependent protein catabolic process', 'transcription from RNA polymerase III promoter', 'nucleotide-excision repair & DNA excision', 'androgen receptor signaling pathway', 'negative regulation of ERBB signaling pathway & negative regulation of epidermal growth factor receptor signaling pathway', 'mRNA modification', 'ubiquitin-dependent SMAD protein catabolic process', 'telomere maintenance via telomerase', 'translational termination', 'transcription initiation from RNA polymerase I promoter', 'DNA deamination', 'very-low-density lipoprotein particle remodeling', 'positive regulation of translational initiation', 'negative regulation of GTPase activity', 'fructose 6-phosphate metabolic process', 'regulation of calcidiol 1-monooxygenase activity']
jaccard_mat=m.plot_jaccard_comp(C,category_mat,G,T,20,kind=1,dict_G=dict_G,SAVE=False,PLOT=False,algorithm='GenGO' if m.kind==1 else 'Bayesian')

n=len(total_annotations)
order_found=range(1,n+1)

#Reorder basic data by size to avoid big circles to cover smaller ones
ind=sorted(range(n), reverse=True, key=lambda k: total_annotations[k])
total_annotations=np.array([total_annotations[ind[i]] for i in xrange(n)])

for j in xrange(len(keys_genes)):
    perturbed_mat[j]=np.array([perturbed_mat[j][ind[i]] for i in xrange(n)])
    mean_gene_level_mat[j]=np.array([mean_gene_level_mat[j][ind[i]] for i in xrange(n)])
    jaccard_mat[j]=np.array([jaccard_mat[j][ind[i]] for i in xrange(n)])

order_found=np.array([order_found[ind[i]] for i in xrange(n)])
pvalues=np.array([pvalues[ind[i]] for i in xrange(n)])
term_names=np.array([term_names[ind[i]] for i in xrange(n)])
C_ind=[C[ind[i]] for i in xrange(n)]


jaccard_comp=jaccard_mat[geneslabel_ind_y][:]
mean_gene_level=mean_gene_level_mat[geneslabel_ind_y][:]
perturbed=perturbed_mat[geneslabel_ind_y][:]

ratio_pert_total=np.array([y*1./x for x,y in zip(total_annotations,perturbed)])

########Basic data created

keys_mat=['Nr annotated and part of X','Mean gene expression (X)','Jaccard Similarity (X)']
values_mat=[perturbed_mat,mean_gene_level_mat,jaccard_mat]
dict_mat={k:v for k,v in zip(keys_mat,values_mat)}

keys=['Nr annotated','Nr annotated and part of X','Order Found','log(P-value)','Percentage perturbed','Mean gene expression (X)','Jaccard Similarity (X)']
values=[total_annotations,perturbed,order_found,pvalues,ratio_pert_total,mean_gene_level,jaccard_comp]
dict_choice={k:v for k,v in zip(keys,values)}

xs=np.array(values[xlabel_ind][:])
ys=np.array(values[ylabel_ind][:])

def x_data(label, UPDATE=True):
    global xs, xlabel_ind
    which_genes(keys_genes[geneslabel_ind_x],True,CALL_DATA_FCT=False)
    xs = dict_choice[label]
    xlabel_ind = keys.index(label)
    if UPDATE:
        update()
def y_data(label):
    global ys, ylabel_ind
    which_genes(keys_genes[geneslabel_ind_y],False,CALL_DATA_FCT=False)
    ys = dict_choice[label.replace('Y','X')]
    ylabel_ind = keys.index(label.replace('Y','X'))
    #print min(ys)-0.05*diff,max(ys)+0.05*diff,diff,max(xs),min(xs)
    update()
def update():
    global browser, browser2, line, man, data
    ax.cla()
    log_categories=[0,1]
    if xlabel_ind not in log_categories:
        diffx=max(xs)-min(xs)
        ax.set_xlim(min(xs)-0.05*diffx,max(xs)+0.05*diffx)
    if ylabel_ind not in log_categories:
        diffy=max(ys)-min(ys)
        ax.set_ylim(min(ys)-0.05*diffy,max(ys)+0.05*diffy)
    ax.set_title(plot_title)
    try:
        if xlabel_ind in log_categories and ylabel_ind in log_categories:
            line, = ax.loglog(xs, ys, 'o', picker=accuracy, visible=False)
        elif xlabel_ind in log_categories:
            line, = ax.semilogx(xs, ys, 'o', picker=accuracy, visible=False)
        elif ylabel_ind in log_categories:
            line, = ax.semilogy(xs, ys, 'o', picker=accuracy, visible=False)
        else:
            line, = ax.plot(xs, ys, 'o', picker=accuracy, visible=False)  # 5 points tolerance
    except ValueError: #sometimes sth weird doesn't allow log plots, ValueError: Data has no positive values, and therefore can not be log-
        line, = ax.plot(xs, ys, 'o', picker=accuracy, visible=False)  # 5 points tolerance

    ax.set_xlabel(keys[xlabel_ind].replace('X',keys_genes[geneslabel_ind_x][:keys_genes[geneslabel_ind_x].find('(')-1]))
    ax.set_ylabel(keys[ylabel_ind].replace('X',keys_genes[geneslabel_ind_y][:keys_genes[geneslabel_ind_y].find('(')-1]))
    browser = PointBrowser(active=browser.active,lastind=browser.lastind,print_values=browser.print_values,otherBrowser=browser2)
    browser2 = PointBrowser(active=browser2.active,lastind=browser2.lastind,indent=0.5,indent2=0.35,color='cyan',print_values=browser2.print_values,otherBrowser=browser)
    browser.otherBrowser = browser2
    #browser2.otherBrowser = browser
    data = [Datum(x,y,z) for x,y,z in zip(xs,ys,total_annotations)]
    man=Manager(ax,ax2,data)
    put_jaccard_indices()
    
    fig.canvas.mpl_connect('pick_event', browser.onpick)
    fig.canvas.mpl_connect('key_press_event', browser.onpress)
    draw()


def command_x(label):
    which_genes(label,1,CALL_DATA_FCT=True)

def command_y(label):
    which_genes(label,0,CALL_DATA_FCT=True)

def which_genes(label,x1_or_y0,CALL_DATA_FCT=False):
    global perturbed, mean_gene_level, jaccard_comp, geneslabel_ind_x, geneslabel_ind_y, values, dict_choice
    if x1_or_y0==1:
        geneslabel_ind_x = keys_genes.index(label)
        geneslabel_ind = geneslabel_ind_x
    else:
        geneslabel_ind_y = keys_genes.index(label)
        geneslabel_ind = geneslabel_ind_y
    if m.s==0:
        m.s=len(G)-1
    jaccard_comp=jaccard_mat[geneslabel_ind][:]
    mean_gene_level=mean_gene_level_mat[geneslabel_ind][:]
    perturbed=perturbed_mat[geneslabel_ind][:]
    
    values=[total_annotations,perturbed,order_found,pvalues,ratio_pert_total,mean_gene_level,jaccard_comp]
    dict_choice={k:v for k,v in zip(keys,values)}
    if CALL_DATA_FCT==True:
        if x1_or_y0==1:
            x_data(keys[xlabel_ind])
        else:
            y_data(keys[ylabel_ind])

def put_jaccard_indices():
    global jaccard_text
    print "Haaaalllo"
    #Jaccard comparison, doesn't update yet when geneslabel are changed
    try:
        choice=[C_ind[i] for i in [browser.lastind,browser2.lastind]]
        print choice,category_mat[geneslabel_ind_x],category_mat[geneslabel_ind_y]
        comp_print=('Jaccard(X): '+str(round(m.jaccard_comp(choice,category_mat[geneslabel_ind_x],G,T,kind=1)[0],5))+'\n'
                    'Jaccard(Y): '+str(round(m.jaccard_comp(choice,category_mat[geneslabel_ind_y],G,T,kind=1)[0],5)))
        jaccard_text.set_text(comp_print)
        fig.canvas.draw()
    except:
        print "But not here",browser.lastind,browser2.lastind
        pass
    

fig = figure(figsize=(14,8))
gs = gridspec.GridSpec(2,1,height_ratios=[5,3])
ax = fig.add_subplot(gs[0])
ax.set_title(plot_title)
#line, = ax.loglog(xs, ys, 'o', picker=accuracy, visible=False)  # 5 points tolerance
ax.set_xlabel(keys[xlabel_ind])
ax.set_ylabel(keys[ylabel_ind])
ax2 = fig.add_subplot(gs[1],xTicks=[],yTicks=[])
subplots_adjust(left=0.35)

axcolor = 'lightgoldenrodyellow'
rax = axes([0.02, 0.7, 0.26, 0.25], title="x-data", axisbg=axcolor)
radio_x = RadioButtons(rax, tuple(keys),active=xlabel_ind)
radio_x.on_clicked(x_data)
rax = axes([0.02, 0.4, 0.26, 0.25], title="y-data", axisbg=axcolor)
radio_y = RadioButtons(rax, tuple([k.replace('X','Y') for k in keys]),active=ylabel_ind)
radio_y.on_clicked(y_data)
rax = axes([0.02, 0.21, 0.26, 0.14], title="which genes (x-axis), X", axisbg=axcolor)
radio_genes_x = RadioButtons(rax, tuple(keys_genes),active=geneslabel_ind_x)
radio_genes_x.on_clicked(command_x)
rax = axes([0.02, 0.02, 0.26, 0.14], title="which genes (y-axis), Y", axisbg=axcolor)
radio_genes_y = RadioButtons(rax, tuple(keys_genes),active=geneslabel_ind_y)
radio_genes_y.on_clicked(command_y)

jaccard_text = ax2.text(0.7, 0.6,  '',transform=ax2.transAxes, va='top')

class Datum:
    colorin = colorConverter.to_rgba('red')
    colorout = colorConverter.to_rgba('white')
    def __init__(self, x, y, size, focus=False):
        self.x = x
        self.y = y
        self.size = size
        if focus: self.color = Datum.colorin
        else: self.color = Datum.colorout

class Manager:
    def __init__(self, ax, ax2, data):
        self.axes = ax
        self.canvas = ax.figure.canvas
        self.data = data
        self.Nxy = len(data)
        self.xys = [(d.x, d.y) for d in data]
        sizes=[d.size for d in data]
        fig = ax.figure
        self.collection = RegularPolyCollection(
            fig.dpi, 6, sizes=tuple(sizes),
            offsets = self.xys,
            facecolors = [d.color for d in data],
            transOffset = ax.transData)
        ax.add_collection(self.collection)
        self.text = ax2.text(0.7, 0.6,  '',
                            transform=ax2.transAxes, va='top')
        

class PointBrowser:
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and pervious points
    """
    def __init__(self,active=True,lastind=None,indent=0.325,indent2=0.5,color='red',print_values='',otherBrowser=None):
        self.lastind = 0
        self.active = active
        self.otherBrowser = otherBrowser
        self.indent = indent
        self.color=color
        self.print_values=print_values
        self.text = ax.text(0.02, indent2,  '', color=self.color,
                            transform=ax.transAxes, va='top')

        if lastind!=None:
            self.lastind=lastind
            self.selected,  = ax.plot([xs[lastind]], [ys[lastind]], 'o', alpha=0.2,
                                      color=self.color)
            self.update(EVERYTHING=False)
        else:
            self.lastind=None
            self.selected,  = ax.plot([xs[0]], [ys[0]], 'o', alpha=0.2,
                                      color=self.color, visible=False)

    def change_activity(self):
        self.active = not self.active
        self.otherBrowser.active = not self.otherBrowser.active

    def onpress(self, event):
        if self.lastind is None: return
        if event.key not in ('n', 'p'): return
        if self.active == False:
            #print self.indent, self.active
            return self.otherBrowser.onpress(event)
        if event.key=='n': inc = 1
        else:  inc = -1


        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(xs)-1)
        self.update()

    def onpick(self, event):
        if event.artist!=line: return True

        N = len(event.ind)
        if not N: return True

        if self.active == False:
            #print self.indent, self.active
            return self.otherBrowser.onpick(event)

        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        try:
            distances = np.hypot(x-xs[event.ind], y-ys[event.ind])

            indmin = distances.argmin()
            if self.lastind in event.ind:
                distances = distances.tolist()
                e=event.ind.tolist()
                min_d=min(distances)
                try:
                    indmin=distances.index(min_d,e.index(self.lastind)+1)
                except ValueError:
                    pass
            dataind = event.ind[indmin]

            self.lastind = dataind
            self.update()
        except TypeError: #happens sometimes at the boundary of ax
            pass

    def update(self, EVERYTHING=True):
        global jaccard_text
        if self.lastind is None: return

        #print self.indent #as identifier

        dataind = self.lastind

        print "Here",self.indent

        if EVERYTHING==True:
            print "But here!"
            ax2.cla()

            to_print=''
            self.print_values=''
            for v,k in zip(values,keys):
                if 'X' in k:
                    helper=dict_mat[k][geneslabel_ind_x][dataind]
                    if type(helper)==type(np.array([1.])[0]):
                        helper=round(helper,3)
                    to_print+=k+': '+'\n'
                    self.print_values+=str(helper)+'\n'
                    helper=dict_mat[k][geneslabel_ind_y][dataind]
                    if type(helper)==type(np.array([1.])[0]):
                        helper=round(helper,3)
                    to_print+=k.replace('X','Y')+': '+'\n'
                    self.print_values+=str(helper)+'\n'
                else:
                    helper=v[dataind]
                    if type(helper)==type(np.array([1.])[0]):
                        helper=round(helper,3)
                    to_print+=k+': '+'\n'
                    self.print_values+=str(helper)+'\n'
            #print to_print

            ax2.text(0.02, 0.97, to_print, va='top')
            ax2.text(self.indent, 0.968, self.print_values, color=self.color,va='top')
            ax2.text(self.otherBrowser.indent, 0.968, self.otherBrowser.print_values, color=self.otherBrowser.color,va='top')

        
        self.selected.set_markersize(math.sqrt(total_annotations[dataind])+10)
        self.selected.set_visible(True)
        self.selected.set_data(xs[dataind], ys[dataind])
        if len(term_names[dataind])>85:
            ind=term_names[dataind][:85].rfind(' ')
            self.text.set_text(term_names[dataind][:ind]+"\n"+term_names[dataind][ind+1:])
        else:
            self.text.set_text(term_names[dataind])
        ax2.set_xticks([])
        ax2.set_yticks([])
        #set_browser_order()
        if EVERYTHING==True:
            print "a"
            self.change_activity()
            put_jaccard_indices()
        fig.canvas.draw()

browser = PointBrowser()
browser2 = PointBrowser(active=False,lastind=None,indent=0.5,indent2=0.35,color='cyan',otherBrowser=browser)
browser.otherBrowser = browser2
#browser2.otherBrowser = browser
man=None
data=None
line=None

x_data(keys[xlabel_ind],False)
y_data(keys[ylabel_ind])

show()
