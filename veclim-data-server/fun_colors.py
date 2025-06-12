import numpy
import matplotlib as mpl
from matplotlib import pyplot as plt

cmaps = {}


clscl = ['#50c0ad','#8dcbc1','#c6e0ee','white','#f5d9b8','#e2988a','#f15a48']
clbins = [-3,-2,-1,1,2,3]
cllbl = ["1/16","1/8","1/4","1/2","2","4","8","16"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='both')
tran = lambda x: x

cmaps['diff_colegg'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}


clscl = ['#50c0ad','#8dcbc1','#c6e0ee','white','#f5d9b8','#e2988a','#f15a48']
clbins = [-0.15,-0.1,-0.05,0.05,0.1,0.15]
cllbl = ['-20%','-15%','-10%','-5%','5%','10%','15%','20%']

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='both')
tran = lambda x: numpy.nanmean(x,axis=2)

cmaps['diff_pouts'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}


clscl = ['#50c0ad','#8dcbc1','#c6e0ee','white','#f5d9b8','#e2988a','#f15a48']
clbins = [-20.0,-10.0,-1.0,1.0,10.0,20.0]
cllbl = ['-30','-20','-10','-1','1','10','20','30']

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='both')
tran = lambda x: numpy.nanmean(x*4000,axis=2)

cmaps['diff_iouts'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}


clscl = ["#f7fcf5", "#d5efcf", "#9ed898", "#54b567", "#1d8641", "#00441b"]
clbins = [0,1,2,3,4,5,6]
cllbl = ["0","10","100","200","300","400","500"]

cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=False, extend='neither')
tran = lambda x: x

cmaps['papatasi'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}


clscl = ['#00000000', '#fbe590', '#f7a034', '#e85229', '#931b1f']
clbins = [0,0.01,0.05,0.1,0.2,0.5]
cllbl = ["0","1","5","10","20","50"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')
tran = lambda x: numpy.nanmean(x,axis=2)

cmaps['pouts'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}


clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#e85229', '#931b1f']
clbins = [0,1,10,50,100,500,1000]
cllbl = ["0","1","10","50","100","500","1000"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')
tran = lambda x: numpy.nanmean(x*4000,axis=2)

cmaps['iouts'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = ["black"]
clbins = [0,1]
cllbl = ["0"]

cmap = None
norm = None
tran = lambda x: x

cmaps['vabun'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = numpy.array([
         '#081d58',
         '#253494',
         '#225ea8',
         '#1d91c0',
         '#41b6c4',
         '#7fcdbb',
         '#fcc65a', 
         '#f7a034', 
         '#f47b2c', 
         '#e85229', 
         '#d82929',
         '#b42125',
         '#00000000',
         ])[::-1].tolist()
clbins = numpy.cumsum([0,1,30,28,31,30,31,30,31,31,30,31,30,31])
cllbl = ["NA","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')
def trann2(mat):
    tmp = mat[:,:,:]>1.0
    annmat = numpy.sum(numpy.cumsum(tmp[:,:,1:]>tmp[:,:,:-1],axis=2)==0,axis=2)+1
    annmat[annmat==365] = 0
    return annmat

cmaps['larva'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': trann2
}

clscl = ['#00000000', '#931b1f']
clbins = [0,0.5,1]
cllbl = ["Unknown/absent","","Reported/established"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')
tran = lambda x: numpy.array(x)

cmaps['presence'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = ['#00000000', '#f15a48', '#1b3958', '#167997', '#50c0ad']
clbins = [0,1,2,3,4,5]
cllbl = ["Unknown/absent","","Global presence (2024)", "VectAbundance (2010-2022)", "AIMsurv (2020)", "VectorBase (2010-2024)"]
cmap = None
norm = None
tran = lambda x: numpy.array(x)

cmaps['albosurv'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#f47b2c', '#e85229', '#d82929', '#931b1f']
clbins = [-4,-3,-2,-1,0,1,2,3,4]
cllbl = ["1/16","1/8","1/4","1/2","1","2","4","8","16"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')
tran = lambda x: numpy.log2(numpy.nanmean(x,axis=2)/5.0) # for '2010-2019' (corrected)

cmaps['FuzzyLocV6'] = {
    'clscl': clscl,
    'clbins': clbins,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = ['#00000000', '#fbe590', '#fcc65a', '#f7a034', '#f47b2c', '#e85229', '#d82929', '#931b1f']
clbins = [-4,-3,-2,-1,0,1,2,3,4]
cllbl = ["NA","Trace","1/2","1","2","4","8","16","32"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')
tran = lambda x: numpy.log2(numpy.nanmean(x,axis=2)/2.1205509) # for '2010-2019'

cmaps['FuzzyLoc'] = {
    'clscl': clscl,
    'clbins': cllbl,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = ['#FFFFFF','#C2C200','#FFE53B','#FFA500','#FF0000','#EE82EE','#A020F0','#4C00FF','#00008B']
clbins = [-4,-3,-2,-1,0,1,2,3,4]
cllbl = ["NA","Trace","1/2","1","2","4","8","16","32"]
cmap = mpl.colors.ListedColormap([mpl.colors.to_rgba(c) for c in clscl])
norm = mpl.colors.BoundaryNorm(clbins, cmap.N, clip=True, extend='neither')
tran = lambda x: numpy.log2(numpy.nanmean(x,axis=2)/0.384307) # for 'Q4.a100+1'

cmaps['PopDM'] = {
    'clscl': clscl,
    'clbins': cllbl,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = []
clbins = []
cllbl = []
cmap = plt.cm.viridis
norm = mpl.colors.Normalize(vmin=-3,vmax=3)
tran = lambda x: numpy.log10(numpy.nanmean(x,axis=2))

cmaps['viridis'] = {
    'clscl': clscl,
    'clbins': cllbl,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}

clscl = []
clbins = []
cllbl = []
cmap = plt.cm.viridis
norm = mpl.colors.Normalize(vmin=-3,vmax=3)
tran = lambda x: numpy.log10(x)

cmaps['viridis.pop'] = {
    'clscl': clscl,
    'clbins': cllbl,
    'cllbl': cllbl,
    'cmap': cmap,
    'norm': norm,
    'tran': tran
}
