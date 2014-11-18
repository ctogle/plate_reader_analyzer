import modular_core.libfundamental as lfu
import modular_core.libpostprocess as lpp
import modular_core.libdatacontrol as ldc

import types
import numpy as np
from copy import deepcopy as dcopy

import pdb

class post_process_plate_reader(lpp.post_process):

    def __init__(self, *args, **kwargs):
        self.impose_default('_always_sourceable_',[],**kwargs)
        lpp.post_process.__init__(self, *args, **kwargs)

    def provide_axes_manager_input(self, lp = True, cp = False, 
            bp = False, vp = False, tp = True, x_title = 'x-title', 
            y_title = 'y-title', title = 'title', plt_callbacks = None):
        self.use_line_plot = lp
        self.use_color_plot = cp
        self.use_bar_plot = bp
        self.use_voxel_plot = vp
        self.use_table_plot = tp
        self.x_title = x_title
        #meas = self._measurement_
        self.y_title = y_title
        self.title = y_title
        if plt_callbacks is None:plt_callbacks = {}
        self.plt_callbacks = plt_callbacks
        #self.plt_callbacks = {
        #    'change_y_domain':self.hello}

    #def hello(self): print 'hey there!'

    def start_pool(self, *args, **kwargs):
        #pool = []
        #if 'simulation' in self.input_regime:
        pool = args[1]
        return pool

    def postproc(self, *args, **kwargs):
        #pool should invariably be a list of 
        #  lists of data objects for each trajectory
        #a method must always be provided by a superclass
        #  a pool and a p_space are optional, default
        #  is to use the ensemble
        self.determine_regime(args[0])
        #pool = []
        pool = self.start_pool(*args, **kwargs)
        sources = self.get_source_reference(1, *args, **kwargs)
        from_always = self.parent.parent._children_[0].data
        if type(self.input_regime) is types.ListType:
            for inp in self.input_regime:
                inpmobj = lfu.grab_mobj_by_name(inp, from_always)
                lfu.zip_list(pool, [inpmobj.data.data])
        else:
            inpmobj = lfu.grab_mobj_by_name(self.input_regime, from_always)
            lfu.zip_list(pool, [inpmobj.data.data])
        for src in sources: lfu.zip_list(pool, src.data)
        if self.regime == 'all trajectories':
            self.handle_all_trajectories(kwargs['method'], pool)
        elif self.regime == 'manual grouping':
            self.handle_manual_grouping(kwargs['method'], pool)
        elif self.regime == 'per trajectory':
            self.handle_per_trajectory(kwargs['method'], pool)

class post_process_extract_data(post_process_plate_reader):

    #def provide_axes_manager_input(self, *args, **kwargs):
    #    cbs = {
    #        'change_y_domain':self.update}
    #    kwargs['plt_callbacks'] = cbs
    #    post_process_plate_reader.provide_axes_manager_input(self,*args,**kwargs)

    def __init__(self, *args, **kwargs):
        if not 'base_class' in kwargs.keys():
            kwargs['base_class'] =\
                lfu.interface_template_class(object, 'extract')
        if not 'label' in kwargs.keys():
            kwargs['label'] = 'ratios'
        if not 'valid_regimes' in kwargs.keys():
            kwargs['valid_regimes'] = ['all trajectories']
        if not 'regime' in kwargs.keys():
            kwargs['regime'] = 'all trajectories'
        kwargs['_single_input_'] = True
        self.impose_default('x_domain', None, **kwargs)
        self.impose_default('rat_domain', None, **kwargs)
        #self.impose_default('function_of', None, **kwargs)
        post_process_plate_reader.__init__(self, *args, **kwargs)

    #def update(self, *args, **kwargs):
    #    combo = args[0]
    #    comdx = combo.currentIndex()
    #    self.rat_domain = combo.itemText(comdx)
    #    updated = self.ratios(*self.raw, **kwargs)
    #    for da in self.data[0]:
    #        up = lfu.grab_mobj_by_name(da.label, updated)
    #        da.scalars = up.scalars

    def postproc(self, *args, **kwargs):
        kwargs['method'] = self.extract
        post_process_plate_reader.postproc(self, *args, **kwargs)

    def extract(self, *args, **kwargs): 
        pool = args[0][0]
        #self.raw = args
        newpool = []
        stdpool = []
        for da in pool:
            if hasattr(da,'subscalars'):
                stdlabel = da.label + '-stddev'
                stddat = da.subscalars[0] - da.scalars
                std = ldc.scalars(label = stdlabel,scalars = stddat)
                stdpool.append(std)
            newpool.append(dcopy(da))
        newpool.extend(stdpool)
        return newpool

    def set_target_settables(self, *args, **kwargs):
        self.valid_regimes = ['all trajectories']
        self.valid_inputs = self.get_valid_inputs(*args, **kwargs)
        always = args[1]._children_[0].data
        always = [al for al in always 
            if al.label == self.input_regime]
            #if al.label in self._always_sourceable_]
        kwargs['_always_sources_'] = always
        capture_targetable = self.get_targetables(*args, **kwargs)
        self.target_list = capture_targetable[:]
        self.capture_targets = self.target_list
        post_process_plate_reader.set_target_settables(self, *args, **kwargs)

    def set_settables(self, *args, **kwargs):
        self.handle_widget_inheritance(*args, from_sub = False)
        #capture_targetable = self.get_targetables(*args, **kwargs)
        exempt = ['Time','T Deg']
        accounted = []
        new = []
        for tg in self.target_list:
            if tg.strip() not in exempt and not tg in accounted:
                #tdx = self.target_list.index(tg)
                #self.target_list.insert(tdx+1,tg+'-stddev')
                new.append(tg+'-stddev')
                accounted.append(tg)
        self.target_list.extend(new)
        print 'new tg list!'
        print self.target_list
        capture_targetable = self.target_list
        post_process_plate_reader.set_settables(
            self, *args, from_sub = True)

class post_process_ratios(post_process_plate_reader):

    def provide_axes_manager_input(self, *args, **kwargs):
        cbs = {
            'change_y_domain':self.update}
        kwargs['plt_callbacks'] = cbs
        post_process_plate_reader.provide_axes_manager_input(self,*args,**kwargs)

    def __init__(self, *args, **kwargs):
        if not 'base_class' in kwargs.keys():
            kwargs['base_class'] =\
                lfu.interface_template_class(object, 'ratios')
        if not 'label' in kwargs.keys():
            kwargs['label'] = 'ratios'
        if not 'valid_regimes' in kwargs.keys():
            kwargs['valid_regimes'] = ['all trajectories']
        if not 'regime' in kwargs.keys():
            kwargs['regime'] = 'all trajectories'
        kwargs['_single_input_'] = True
        self.impose_default('x_domain', None, **kwargs)
        self.impose_default('rat_domain', None, **kwargs)
        #self.impose_default('function_of', None, **kwargs)
        post_process_plate_reader.__init__(self, *args, **kwargs)

    def update(self, *args, **kwargs):
        combo = args[0]
        comdx = combo.currentIndex()
        self.rat_domain = combo.itemText(comdx)
        updated = self.ratios(*self.raw, **kwargs)
        for da in self.data[0]:
            up = lfu.grab_mobj_by_name(da.label, updated)
            da.scalars = up.scalars

    def postproc(self, *args, **kwargs):
        kwargs['method'] = self.ratios
        post_process_plate_reader.postproc(self, *args, **kwargs)

    def ratios(self, *args, **kwargs):
        pool = args[0][0]
        self.raw = args
        trudom = self.x_domain
        rat_domain = lfu.grab_mobj_by_name(self.rat_domain, pool)
        rat_dom_std = rat_domain.subscalars[0] - rat_domain.scalars
        rdomv = rat_domain.scalars
        dlabs = [d.label for d in pool]
        data = ldc.scalars_from_labels(dlabs)
        stddata = []
        def div(u,v):return v/u if not u == 0.0 else float(10000000000)
        for dx in range(len(data)):
            if data[dx].label == trudom:
                data[dx].scalars = pool[dx].scalars[:]
            else:
                data[dx].scalars = np.array([div(u,v) for v,u in 
                    zip(pool[dx].scalars,rdomv)], dtype = float)
                #if hasattr(rat_domain,'subscalars') and hasattr(pool[dx],'subscalars'):
                #    pstd = pool[dx].subscalars[0] - pool[dx].scalars
                #    pdb.set_trace()
                
                # NOT CORRECT stdd = np.std(data[dx].scalars)
                #stddata.append(ldc.scalars(
                #    label = data[dx].label+'-stddev', scalars = stdd))
        
        #pdb.set_trace()
        data.extend(stddata)
        return data

    def set_target_settables(self, *args, **kwargs):
        self.valid_regimes = ['all trajectories']
        self.valid_inputs = self.get_valid_inputs(*args, **kwargs)
        always = args[1]._children_[0].data
        always = [al for al in always 
            if al.label == self.input_regime]
            #if al.label in self._always_sourceable_]
        kwargs['_always_sources_'] = always
        capture_targetable = self.get_targetables(*args, **kwargs)
        self.target_list = capture_targetable[:]
        self.capture_targets = self.target_list
        post_process_plate_reader.set_target_settables(self, *args, **kwargs)

    def set_settables(self, *args, **kwargs):
        self.handle_widget_inheritance(*args, from_sub = False)
        #capture_targetable = self.get_targetables(*args, **kwargs)
        capture_targetable = self.target_list
        if self.x_domain is None:
            xdom = capture_targetable[0]
            self.x_domain = xdom
        else: xdom = self.x_domain
        if self.rat_domain is None:
            ratdom = capture_targetable[0]
            self.rat_domain = ratdom
        else: ratdom = self.rat_domain
        self.widg_templates.append(
            lgm.interface_template_gui(
                widgets = ['radio'], 
                labels = [capture_targetable], 
                instances = [[self]], 
                keys = [['rat_domain']], 
                box_labels = ['Ratio Denominator'], 
                initials = [[ratdom]]))
        self.widg_templates.append(
            lgm.interface_template_gui(
                widgets = ['radio'], 
                labels = [capture_targetable], 
                instances = [[self]], 
                keys = [['x_domain']], 
                box_labels = ['X-Domain (Exempt)'], 
                initials = [[xdom]]))
        #self.x_domain = 'Time'
        #self.rat_domain = pool[2].label
        post_process_plate_reader.set_settables(
            self, *args, from_sub = True)

valid_postproc_base_classes = [
    lfu.interface_template_class(
        post_process_ratios, 'ratios'),
    lfu.interface_template_class(
        post_process_extract_data, 'extract')]
lpp.valid_postproc_base_classes = valid_postproc_base_classes

if __name__ == 'modular_pra.libplatereaderprocesses':
  if lfu.gui_pack is None: lfu.find_gui_pack()
  lgm = lfu.gui_pack.lgm
  lgd = lfu.gui_pack.lgd
  lgb = lfu.gui_pack.lgb

if __name__ == '__main__':
  print 'this is a library!'





