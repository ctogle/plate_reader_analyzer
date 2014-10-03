import modular_core.libfundamental as lfu
import modular_core.libsettings as lset
import modular_core.liboutput as lo
import modular_core.libmath as lm
#import modular_core.libgeometry as lgeo
import modular_core.libdatacontrol as ldc
import modular_core.libpostprocess as lpp
#import libs.plate_reader_analyzer.libplatereaderprocesses as lpdap

from collections import OrderedDict
import numpy as np
from scipy.optimize import curve_fit as cufit

from copy import deepcopy as dcopy
import os, sys, time, traceback, re
import pdb



if __name__ == 'modular_pra.libplate_reader_analyzer':
	if lfu.gui_pack is None: lfu.find_gui_pack()
	lgm = lfu.gui_pack.lgm
	lgd = lfu.gui_pack.lgd
	lgb = lfu.gui_pack.lgb

if __name__ == '__main__': print 'this is a library!'

class data_block(lfu.modular_object_qt):
	swaps = {
            '\xb0':' Deg ',
            '\r':'',
                }
	def __init__(self, *args, **kwargs):
		self.raw = args[0]
		lfu.modular_object_qt.__init__(self, *args, **kwargs)
	def display_read(self):
		for key in self.read.keys():
			print '\nkey :', key
			for li in self.read[key]: print '\t' + li
	def _read_(self, *args, **kwargs):
		self.read = args[0]
		#self.display_read()
	def swap_check(self, li):
		return swap_check(self.swaps, li)

def swap_check(swaps, li):
    for sw in swaps.keys():
        li = li.replace(sw, swaps[sw])
    return li

class header_block(data_block):
	def __init__(self, *args, **kwargs):
		blks = args[0]
		#try:
		#	raw = blks[0].raw + blks[1].raw + blks[2].raw
		#except IndexError: pdb.set_trace()
                raw = []
                for bl in blks: raw += bl.raw#THIS WILL BE STUPID SLOW...
		data_block.__init__(self, raw, **kwargs)
	def _read_(self, pra):
		read = OrderedDict()
		for li in self.raw:
			ke = li[:li.find(',')].replace(':','')
			read[ke] = [li[li.find(',')+1:]]
		data_block._read_(self, read)

class procedure_block(data_block):
	def __init__(self, *args, **kwargs):
		blks = args[0]
		#raw = blks[0].raw + blks[1].raw
                raw = []
                for bl in blks: raw += bl.raw#THIS WILL BE STUPID SLOW...
		data_block.__init__(self, raw, **kwargs)
	def _read_(self, pra):
		read = OrderedDict()
		relev = self.raw[1:]
		kins = []; kin_on = False
		for li_dex in range(len(relev)):
			li = self.swap_check(relev[li_dex])
			sp = [x.strip() for x in li.split(',')]
			ke = sp[0]
			if ke.startswith('End Kinetic'):
				kin_on = False
			elif kin_on: read[kin_ke].append(li.strip())
			elif li.startswith(','):
				read[read.keys()[-1]].append(li.strip())
			elif ke.startswith('Start Kinetic'):
				kins.append(li_dex); kin_on = True
				kin_ke = 'Kinetic #' + str(len(kins))
				read[kin_ke] = [','.join(sp[1:])]
			else: read[ke] = [','.join(sp[1:])]
		data_block._read_(self, read)

class layout_block(data_block):
	def __init__(self, *args, **kwargs):
		blks = args[0]
		raw = args[0].raw
		data_block.__init__(self, raw, **kwargs)

	def _read_(self, pra):
		sps = [r.split(',') for r in self.raw[2:]]
		rows = len(sps)
		cols = len(sps[1])-1
		reduced = np.ndarray((rows, cols), dtype = object)
		for row in range(rows): reduced[row,:] = sps[row][:-1]
		read = OrderedDict()
		read['rows'] = self.raw[1].split(',')[1:]
		read['cols'] = reduced[:,0]
		read['table'] = reduced[:,1:]
		data_block._read_(self, read)

class obs_data_block(data_block):
	_modu_program_ = 'plate_reader_analyzer'

	def __init__(self, *args, **kwargs):
                self.impose_default('is_OD_block', False, **kwargs)
		self.impose_default('capture_targets',[],**kwargs)
		self.impose_default('replicate_reduced',False,**kwargs)
		self.impose_default('normalized_reduced',False,**kwargs)
		self.impose_default('timept_filtered',False,**kwargs)
		self.impose_default('timept_filtered_OD',False,**kwargs)
		self.impose_default('background_subtracted',False,**kwargs)
                self.impose_default('blank_well_filter_std_factor',5,**kwargs)
                self.impose_default('fake_zero_value',0.000000001,**kwargs)
		self.impose_default('phase_reduced', False, **kwargs)
                self.impose_default('phase_type', 'Log', **kwargs) 
                self.impose_default('override_domains_with_OD', False, **kwargs)
                self.impose_default('override_thresholds', False, **kwargs)
		blks = args[0]
		raw = blks[0].raw + blks[1].raw
		data_block.__init__(self, *(raw,), **kwargs)
		if not hasattr(self, 'label'): self.label = kwargs['label']
		self.output = lo.output_plan(
			label = self.label + ' Output', parent = self)
		self.output.flat_data = True
		self._children_ = [self.output]
		self.output.output_plt = False

        def determine_replicates(self, pra, wells):
            repcnt = int(pra.reps_per_rep)
            grouped = [wells[k::repcnt] for k in range(repcnt)]
            grouped = zip(*grouped)
            repkey = []
            for gr in grouped:
                #nums = [re.search('\d+', g) for g in gr]
                repkey.append(';'.join(gr))
            return repkey

	def _read_(self, pra):
		meas = self.swap_check(self.raw[0])
		self._measurement_ = meas
		layout = self.swap_check(self.raw[1]).split(meas)
                try:wells = layout[1].replace(',', '',1).split(',')
                except:pdb.set_trace()
		well_count = len(wells)
		self._well_key_ = wells
                pra._well_key_ = wells
                repkey = self.determine_replicates(pra, wells)
                self._replicate_key_ = repkey
                pra._replicate_key_ = repkey
		conds = layout[0].split(',')
		cond_count = len(conds)
		self._cond_key_ = conds
		data_lines = self.raw[2:]
		meas_count = len(data_lines)
		_cond_data_ = np.ndarray((meas_count,cond_count),dtype=object)
		_well_data_ = np.ndarray((meas_count,well_count),dtype=object)
                self.replicate_mobjs = {}
                self.well_mobjs = {}
                self.cond_mobjs = {}
		for re_dex, re in enumerate(data_lines):
		        sp = re.split(',')
			_cond_data_[re_dex,:] = sp[:cond_count]
			_well_data_[re_dex,:] = sp[cond_count:]
                for codex, cond in enumerate(self._cond_key_):
                        self.cond_mobjs[cond] = cond_data(cond)
                if self.is_OD_block:
                    lo = self.OD_threshold_low
                    mi = self.OD_threshold_middle
                    hi = self.OD_threshold_high
                else: lo,mi,hi = None, None, None
                for wedex, well in enumerate(self._well_key_):
                        self.well_mobjs[well] =\
                            well_data(well, od_data = self.is_OD_block, 
                                lo_thresh = lo, mi_thresh = mi, hi_thresh = hi, 
                                parent = self)
                for repdex, rep in enumerate(self._replicate_key_):
                    repwells = []
                    for wekey in self._well_key_:
                        if rep.count(wekey) > 0:
                            repwells.append(self.well_mobjs[wekey])
                    self.replicate_mobjs[rep] = replicate_data(
                        rep, od_data = self.is_OD_block, 
                        parent = self, wells = repwells)
		self._cond_data_ = _cond_data_
		self._well_data_ = _well_data_
		read = OrderedDict()
		read['measurement'] = [self._measurement_]
		read['condition_key'] = self._cond_key_
		read['condition_data'] = self._cond_data_
		read['wells_key'] = self._well_key_
		read['wells_data'] = self._well_data_
	        read['well_mobjs'] = self.well_mobjs
	        read['cond_mobjs'] = self.cond_mobjs
		data_block._read_(self, read)

	def deep_parse(self, pra):
		def empty_check(unic):
			checked = [x for x in unic if not x == '']
			return checked
		def ovrflw_check(unic):
			checked =\
				['1000000.0' if x.count('OVRFLW') else x for x in unic]
			return checked
                def back_r_swap(unic):
                        checked = [x.replace('\r', '') for x in unic]
                        return checked
		def _all_checks_(unic):
                        unic = back_r_swap(unic)
			unic = empty_check(unic)
			checked = ovrflw_check(unic)
			return checked
		def hms_to_mins(hms):
			convs = [60.0, 1.0, 1.0/60.0]
			mins = [sum([float(pair[0])*pair[1] for pair 
				in zip(hms[x].split(':'), convs)]) 
				for x in range(len(hms))]
			return mins
		def time_filter(hms):
			hms0 = '0:00:00'
			hms = [min_ for d, min_ in 
				enumerate(hms) if not 
				min_ == hms0 or d == 0]
			mins = hms_to_mins(hms)
			return mins
		known_filters = OrderedDict()
		#known_filters['Time'] = time_filter
		known_filters['_all_'] = _all_checks_
		def filter_(key, dat):
			kf = known_filters
			dat = kf['_all_'](dat)
			if key in kf.keys(): dat = kf[key](dat)
			try: return np.array(dat,dtype='f')
			except ValueError: pdb.set_trace()
		measur = self._measurement_
		weldat = self._well_data_
		welkey = self._well_key_
		condat = self._cond_data_
		conkey = self._cond_key_
                all_data = []
		#con_data = ldc.scalars_from_labels(conkey)
		for dex, key in enumerate(conkey):
                        new = ldc.scalars(label = key, 
                            scalars = filter_(key, condat[:,dex]))
                        self.cond_mobjs[key].data = new
                        self.cond_mobjs[key].process()
                        all_data.append(new)
			#con_data[dex].scalars = filter_(key, condat[:,dex])
		#wel_data = ldc.scalars_from_labels(welkey)
                timedata = self.cond_mobjs['Time'].data
		for dex, key in enumerate(welkey):
                        new = ldc.scalars(label = key, 
                            scalars = filter_(key, weldat[:,dex]))
                        self.well_mobjs[key].data = new
                        self.well_mobjs[key].process(timedata)
                        all_data.append(new)
			#wel_data[dex].scalars = filter_(key, weldat[:,dex])
                repkey = self._replicate_key_
                for dex, key in enumerate(repkey):
                    self.replicate_mobjs[key].process()
                self._unreduced_ = lfu.data_container(data = all_data[:])
                self._all_data_ = all_data
                #self.data = self._unreduced_
		self.apply_reductions()
                #self._reduced_ = self.apply_reduction(self._unreduced_.data)
		#self.update_replicate_reduction()


                '''
                blanks = lfu.uniqfy(pra.template_read['Blank Well'])
                if not len(blanks) == 1:
                    print 'blank wells should not change identity!!'
                    pdb.set_trace()
                spl = blanks[0].split(' ')
                rng = lm.make_range(spl[1])[0].split(',')
                blanks = [''.join([spl[0], str(int(float(r)))]) for r in rng]
                '''
                rpr = int(pra.reps_per_rep)
                if pra.first_blank_well is None:
                    pra.first_blank_well = welkey[-1*rpr]
                blanks = welkey[-1*rpr:]
                self._blank_well_key_ = blanks
                self.blank_well_filter(pra, blanks)

        def blank_well_filter(self, pra, blanks):
            # determine if every blank well data point is an outlier
            # if at a given timepoint, one of 3 blank wells is an outlier
            #  do not use this value to calculate background noise
            # if at a given timepoint, two or more of 3 blank wells are outliers
            #  throw out all data associated with that time point
            # also identify what fraction of each blank well is an outlier and
            #  report

            def outlier(val,mean,std):
                stdthresh = self.blank_well_filter_std_factor*std
                out = abs(val - mean) > stdthresh
                return out

            blank_wells = [self.well_mobjs[ke] for ke in blanks]
            self.blank_well_mobjs = blank_wells
            if not blank_wells:
                print 'NO BLANK WELL DATA IDENTIFIED!!'
                pdb.set_trace()
            means = [d.mean for d in blank_wells]
            stddevs = [d.stdv for d in blank_wells]
            dcount = len(blank_wells[0].data.scalars)

            flags = [[None]*len(blank_wells)]*dcount
            for fldx in xrange(dcount):
                pool = [d.data.scalars[fldx] for d in blank_wells]
                flags[fldx] = [outlier(p,mean,std) for 
                    p,mean,std in zip(pool,means,stddevs)]

            timept_flags = []
            bground_values = []
            for fdx, flgs in enumerate(flags):
                bground_values.extend(
                    [d.data.scalars[fdx] for d,fl in 
                        zip(blank_wells,flgs) if not fl])
                if flgs.count(True) >= 2:timept_flags.append(fdx)
            if not bground_values:
                print 'NO BLANK WELLS FOUND ACCEPTABLE!!'
                pdb.set_trace()

            self.blankwell_outliers = {'total':float(len(flags))}
            tot = self.blankwell_outliers['total']
            for blnk,blankwell in zip(blanks, zip(*flags)):
                badfrac = (blankwell.count(True)/tot)*100.0
                self.blankwell_outliers[blnk] = badfrac

            self.bground_noise_mean = np.mean(bground_values)
            self.bground_noise_stddev = np.std(bground_values)
            self.flagged_timepts = timept_flags

        def apply_normalization_RFU(self, unred):
            if self.is_OD_block: return unred
            OD_block = self.parent.get_OD_block()
            od_data = OD_block.data
            for d, odb in zip(unred.data, od_data.data):
                if d.label in self._well_key_:
                    subd = [x/y for x,y in zip(d.scalars,odb.scalars)] 
                    d.scalars = np.array(subd)
                else: pass
                #else: d.scalars = odb.scalars[:]
            return unred

        def apply_versus_od_reduction(self, unred):
            od_block = self.parent.get_OD_block()
            od_data = lfu.data_container(data = dcopy(od_block._all_data_))
            od_data = od_block.reduce_data(od_data,False, 
                self.timept_filtered, self.background_subtracted, 
                self.phase_reduced, self.replicate_reduced, 
                False, False)

            for d,dod in zip(unred.data,od_data.data):
                if d.label in self._well_key_:
                    d.override_domain = True
                    d.domain = dod.scalars
                    if len(d.domain) != len(d.scalars):
                        pdb.set_trace()
            return unred

        def apply_timept_flags_OD(self, unred):
            pra = self.parent.parent
            blanks = [mobj.well_id for mobj in self.blank_well_mobjs]
            self.blank_well_filter(pra, blanks)
            if self.is_OD_block: return unred
            OD_block = self.parent.get_OD_block()
            flagged_timepts = OD_block.flagged_timepts
            for d in unred.data:
                subd = [d.scalars[vdx] for vdx in xrange(len(d.scalars)) 
                        if not vdx in flagged_timepts]
                d.scalars = subd
            return unred

        def apply_timept_flags(self, unred):
            pra = self.parent.parent
            blanks = [mobj.well_id for mobj in self.blank_well_mobjs]
            self.blank_well_filter(pra, blanks)
            for d in unred.data:
                subd = [d.scalars[vdx] for vdx in xrange(len(d.scalars)) 
                        if not vdx in self.flagged_timepts]
                d.scalars = subd
            return unred

        #consider broken!
	def apply_replicate_reduction(self, unred):
		read = self.parent.parent.read['layout'].read
		flat = lfu.flatten(read['table'])
		well_cnt = len(flat)
		reduced = unred.data[:len(unred.data)-well_cnt]	#list of replicate averaged scalers
		con_offset = len(reduced)
		uniq = lfu.uniqfy(flat)
		layout = OrderedDict()
		for dex, key in enumerate(flat):
			if key in layout.keys(): layout[key].append(dex + con_offset)
			else: layout[key] = [dex + con_offset]
		new = ldc.scalars_from_labels(layout.keys())
		for ndex, key in enumerate(layout.keys()):
			rel_dexes = layout[key]
			rel_dater = [unred.data[d] for d in rel_dexes]
			#rel_dater = [unred[d] for d in rel_dexes]
			zi = zip(*[r.scalars for r in rel_dater])
			new[ndex].scalars = np.array([np.mean(z) for z in zi])
		reduced.extend(new)
		red = lfu.data_container(data = reduced)
		return red

        def apply_background_subtraction(self, unred):
            def smart_subtract(val):
                if val < 0.0: return self.fake_zero_value
                return val
            bg = self.bground_noise_mean
            for d in unred.data:
                if d.label in self._well_key_:
                    subd = np.array([smart_subtract(v - bg) for v in d.scalars])
                    d.scalars = subd
            return unred

        def apply_phase_reduction(self, unred):
            od_block = self.parent.get_OD_block()
            lo_thresh = od_block.OD_threshold_low
            mi_thresh = od_block.OD_threshold_middle
            hi_thresh = od_block.OD_threshold_high
            lo_dex,mi_dex,hi_dex =\
                od_block.get_lo_hi_threshold_indexes(
                    lo_thresh, mi_thresh, hi_thresh)
            if self.phase_type == 'Lag':
                start = 0
                stop = lo_dex
            elif self.phase_type == 'Log':
                start = lo_dex
                stop = mi_dex
            elif self.phase_type == 'Stationary':
                start = mi_dex
                stop = hi_dex
            elif self.phase_type == 'Death':
                start = hi_dex
                stop = None
            else:
                print 'UNKNOWN PHASE TYPE CHOICE', self.phase_type
                return unred
            #filtered = []
            for d in unred.data:
                subd = d.scalars[start:stop]
                d.scalars = subd
                #filt = ldc.scalars(label = d.label, scalars = subd)
                #filtered.append(filt)
            return unred
            #return lfu.data_container(data = filtered)

        def reduce_data(self, unred, tf_od_f, tf_f, 
                bgs_f, phr_f, rred_f, nred_f, domo_f):
            unredtot = len(unred.data[0].scalars)/100.0
            if tf_od_f:
                unredlen = len(unred.data[0].scalars)
                unred = self.apply_timept_flags_OD(unred)
                unredlenpost = len(unred.data[0].scalars)
                self.tf_od_f_cutout = (unredlen - unredlenpost) / unredtot
            else: self.tf_od_f_cutout = 0.0
            if tf_f:
                unredlen = len(unred.data[0].scalars)
                unred = self.apply_timept_flags(unred)
                unredlenpost = len(unred.data[0].scalars)
                self.tf_f_cutout = (unredlen - unredlenpost) / unredtot
            else: self.tf_f_cutout = 0.0
            if bgs_f:
                #unredlen = len(unred.data[0].scalars)
                unred = self.apply_background_subtraction(unred)
                #unredlenpost = len(unred.data[0].scalars)
                #self.bgs_cutout = (unredlen - unredlenpost) / unredtot
            #else: self.bgs_cutout = 0.0
            if phr_f:
                unredlen = len(unred.data[0].scalars)
                unred = self.apply_phase_reduction(unred)
                unredlenpost = len(unred.data[0].scalars)
                self.phr_f_cutout = (unredlen - unredlenpost) / unredtot
            else: self.phr_f_cutout = 0.0
            if rred_f:
                #unredlen = len(unred.data[0].scalars)
                unred = self.apply_replicate_reduction(unred)
                #unredlenpost = len(unred.data[0].scalars)
                #self.rred_f_cutout = (unredlen - unredlenpost) / unredtot
            #else: self.rred_f_cutout = 0.0
            if nred_f:
                #unredlen = len(unred.data[0].scalars)
                unred = self.apply_normalization_RFU(unred)
                #unredlenpost = len(unred.data[0].scalars)
                #self.nred_f_cutout = (unredlen - unredlenpost) / unredtot
            #else: self.nred_f_cutout = 0.0
            if domo_f:
                #unredlen = len(unred.data[0].scalars)
                unred = self.apply_versus_od_reduction(unred)
                #unredlenpost = len(unred.data[0].scalars)
                #self.domo_f_cutout = (unredlen - unredlenpost) / unredtot
            #else: self.domo_f_cutout = 0.0
            return unred

        def apply_reductions(self):
            unred = lfu.data_container(data = dcopy(self._all_data_))
            unred = self.reduce_data(unred, self.timept_filtered_OD, 
                self.timept_filtered, self.background_subtracted, 
                self.phase_reduced, self.replicate_reduced, 
                self.normalized_reduced, self.override_domains_with_OD)
            self.data = unred
            self.capture_targets = [x.label for x in self.data.data]
            self.output.rewidget(True)

	def provide_axes_manager_input(self, 
                    lp = True, cp = False, bp = False, vp = False, tp = True, 
                    x_title = 'x-title', y_title = 'y-title', title = 'title'):
		self.use_line_plot = lp
		self.use_color_plot = cp
		self.use_bar_plot = bp
		self.use_voxel_plot = vp
                self.use_table_plot = tp
		self.x_title = x_title
		meas = self._measurement_
		self.y_title = meas
		self.title = meas

	def set_settables(self, *args, **kwargs):
		window = args[0]
		self.handle_widget_inheritance(*args, **kwargs)
		self.widg_templates.append(                     
			lgm.interface_template_gui(
                                layout = 'horizontal', 
				widgets = ['check_set', 'radio'], 
				labels = [['Apply Phase Reduction'], 
                                    ['Lag', 'Log', 'Stationary', 'Death']], 
				append_instead = [False, None], 
                                refresh = [None,[True]], 
                                window = [None,[window]], 
				instances = [[self], [self]], 
                                initials = [None, [self.phase_type]], 
				keys = [['phase_reduced'], ['phase_type']], 
				callbacks = [
                                    [lgb.create_reset_widgets_wrapper(
			                window, self.apply_reductions)], 
                                    [lgb.create_reset_widgets_wrapper(
			                window, self.apply_reductions)]]))
                self.widg_templates[-1] += lgm.interface_template_gui(
                    widgets = ['text'], 
                    read_only = [True], 
                    box_labels = ['% Of Data Removed By Phase Reduction'], 
                    initials = [[self.phr_f_cutout]])
		self.widg_templates.append(                     
			lgm.interface_template_gui(
                                layout = 'horizontal', 
				widgets = ['check_set'], 
				labels = [['Apply Background Subtraction']], 
				append_instead = [False], 
                                #refresh = [[True]], 
                                #window = [[window]], 
				instances = [[self]], 
				keys = [['background_subtracted']], 
				callbacks = [[lgb.create_reset_widgets_wrapper(
			            window, self.apply_reductions)]]))
                self.widg_templates[-1] += lgm.interface_template_gui(
                    widgets = ['text'], 
                    read_only = [True], 
                    box_labels = ['Mean Background Noise'], 
                    initials = [[self.bground_noise_mean]])
                self.widg_templates[-1] += lgm.interface_template_gui(
                    widgets = ['text'], 
                    read_only = [True], 
                    box_labels = ['Stddev Background Noise'], 
                    initials = [[self.bground_noise_stddev]])
                if not self.is_OD_block:
                    self.widg_templates.append(                     
			lgm.interface_template_gui(
                                layout = 'horizontal', 
				widgets = ['check_set'], 
				labels = [['Apply OD Outlier Reduction']], 
				append_instead = [False], 
                                #refresh = [[True]], 
                                #window = [[window]], 
				instances = [[self]], 
				keys = [['timept_filtered_OD']], 
				callbacks = [[lgb.create_reset_widgets_wrapper(
			            window, self.apply_reductions)]]))
                    self.widg_templates[-1] += lgm.interface_template_gui(
                            widgets = ['text'], 
                            read_only = [True], 
                            box_labels = ['% Of Data Removed By Inheriting OD Outlier Removal'], 
                            initials = [[self.tf_od_f_cutout]])
		    self.widg_templates.append(
			lgm.interface_template_gui(
                                layout = 'horizontal', 
				widgets = ['check_set'], 
				labels = [['Using OD As X-Domain']], 
				append_instead = [False], 
                                #refresh = [[True]], 
                                #window = [[window]], 
				instances = [[self]], 
				keys = [['override_domains_with_OD']], 
				callbacks = [[lgb.create_reset_widgets_wrapper(
				    window, self.apply_reductions)]]))
                    #self.widg_templates[-1] += lgm.interface_template_gui(
                    #        widgets = ['text'], 
                    #        read_only = [True], 
                    #        initials = [[self.domo_f_cutout]])
		    self.widg_templates.append(
			lgm.interface_template_gui(
                                layout = 'horizontal', 
				widgets = ['check_set'], 
				labels = [['Apply Normalization Reduction']], 
				append_instead = [False], 
                                #refresh = [[True]], 
                                #window = [[window]], 
				instances = [[self]], 
				keys = [['normalized_reduced']], 
				callbacks = [[lgb.create_reset_widgets_wrapper(
				    window, self.apply_reductions)]]))
                    #self.widg_templates[-1] += lgm.interface_template_gui(
                    #        widgets = ['text'], 
                    #        read_only = [True], 
                    #        initials = [[self.nred_f_cutout]])
		self.widg_templates.append(                     
			lgm.interface_template_gui(
                                layout = 'horizontal', 
				widgets = ['check_set'], 
				labels = [['Apply Outlier Reduction']], 
				append_instead = [False], 
                                #refresh = [[True]], 
                                #window = [[window]], 
				instances = [[self]], 
				keys = [['timept_filtered']], 
				callbacks = [[lgb.create_reset_widgets_wrapper(
			            window, self.apply_reductions)]]))
                self.widg_templates[-1] += lgm.interface_template_gui(
                        widgets = ['spin'], 
                        keys = [['blank_well_filter_std_factor']], 
                        instances = [[self]], 
                        box_labels = ['Sigma Threshold For Outlier Removal'], 
                        doubles = [[True]], 
                        single_steps = [[1.0]], 
                        initials = [[self.blank_well_filter_std_factor]])
                blnk_well_ids = [we.well_id for we in self.blank_well_mobjs]
                badfraclabels = [
                        '% Of Data Identified As Outliers In Well ' + well 
                        for well in blnk_well_ids]
                badfracvalues = [[self.blankwell_outliers[ke]] 
                                    for ke in blnk_well_ids]
                self.widg_templates[-1] += lgm.interface_template_gui(
                        widgets = ['text','text','text','text'], 
                        read_only = [True, True, True, True], 
                        box_labels = ['% Of Data Removed By Outlier Removal'] +\
                                        badfraclabels, 
                        initials = [[self.tf_f_cutout]]+[b for b in badfracvalues])
		self.widg_templates.append(
			lgm.interface_template_gui(
                                layout = 'horizontal', 
				widgets = ['check_set'], 
				labels = [['Apply Replicate Reduction']], 
				append_instead = [False], 
                                #refresh = [[True]], 
                                #window = [[window]], 
				instances = [[self]], 
				keys = [['replicate_reduced']], 
				callbacks = [[lgb.create_reset_widgets_wrapper(
				    window, self.apply_reductions)]]))
                #self.widg_templates[-1] += lgm.interface_template_gui(
                #        widgets = ['text'], 
                #        read_only = [True], 
                #        initials = [[self.rred_f_cutout]])
		lfu.modular_object_qt.set_settables(
			self, *args, from_sub = True)

class cond_data(object):
	def __init__(self, *args, **kwargs):
		self.cond_id = args[0]

        def process(self):
            self.mean = np.mean(self.data.scalars)
            self.stdv = np.std(self.data.scalars)

def expo(x,a,b,c):
    exp = a * np.exp(b * x) + c
    return exp

def sigmoid(x, a, k, q, b, m, nu):
    if not nu > 0:
        return np.zeros(len(x))
    sig = a + ((k - a)/(1 + q*np.exp(-b*(x - m)))**(1/nu))
    return sig

class replicate_data(object):

    def __init__(self, *args, **kwargs):
        self.replicate_id = args[0]
        self.parent = kwargs['parent']
        self.wells = kwargs['wells']
        try: self.significant_figures = kwargs['significant_figures']
        except KeyError: self.significant_figures = 4
        self.is_OD_data = kwargs['od_data']

    def process(self):
        sfigs = self.significant_figures
        if self.is_OD_data:
            dtimevals = []
            gratevals = []
            for we in self.wells:
                dtimevals.append(we.dtime)
                gratevals.append(we.grate)
            self.mean_doubling_time = np.round(np.mean(dtimevals), sfigs)
            self.mean_growth_rate = np.round(np.mean(gratevals), sfigs)
            self.stddev_doubling_time = np.round(np.std(dtimevals), sfigs)
            self.stddev_growth_rate = np.round(np.std(gratevals), sfigs)

class well_data(object):

	def __init__(self, *args, **kwargs):
            self.parent = kwargs['parent']
	    self.well_id = args[0]
            try: self.significant_figures = kwargs['significant_figures']
            except KeyError: self.significant_figures = 4
            self.is_OD_data = kwargs['od_data']
            if self.is_OD_data:
                self.lo_thresh = kwargs['lo_thresh']
                self.mi_thresh = kwargs['mi_thresh']
                self.hi_thresh = kwargs['hi_thresh']

        def process(self, t):
            self.mean = np.mean(self.data.scalars)
            self.stdv = np.std(self.data.scalars)
            if self.is_OD_data:
                self.dtime, self.grate = self.doubling_time(t.scalars)
                print 'dtimegrate', self.dtime, self.grate

        # find the dexes in the sigmoid fit! its monotonic, smooth, 1-1
        def find_threshold_indices(self,lo,mi,hi,vals):
	    lodelts = [abs(x-lo) for x in vals]
	    lodex = lodelts.index(min(lodelts))
	    midelts = [abs(x-mi) for x in vals]
	    midex = midelts.index(min(midelts))
	    hidelts = [abs(x-hi) for x in vals]
	    hidex = hidelts.index(min(hidelts))
            return lodex,midex,hidex

	def doubling_time(self, t):
            vals = self.data.scalars
            lo = self.lo_thresh
            mi = self.mi_thresh
            hi = self.hi_thresh
            if self.parent.override_thresholds:
                bg = self.parent.bground_noise_mean +\
                    self.parent.bground_noise_stddev
                lo = bg
                hi = max(vals)

	    siglodelts = [abs(x-lo) for x in vals]
	    siglodex = siglodelts.index(min(siglodelts))
            sigdelts = [abs(x-hi) for x in vals]
            sighidex = sigdelts.index(min(sigdelts))
            sigxrelev = t[siglodex:sighidex]/60.0
            sigyrelev = vals[siglodex:sighidex]
            try:
                poptsi, pcovsi = cufit(sigmoid, sigxrelev, sigyrelev)
                sigdat = sigmoid(sigxrelev, *poptsi)
                sigdatdom = 60.0*sigxrelev
                self.cufit_data_sig = ldc.scalars(label = 'sigmoid-fit', 
                    scalars = sigdat, domain = sigdatdom, color = 'r',  
                    override_domain = True)
                if self.parent.override_thresholds:
                    sigmoid_2nd_deriv = [
                        abs(lm.calc_2nd_deriv(sigdatdom, sigdat, k)) 
                            for k in range(1,len(sigdatdom)-1)]
                    infl_deriv = min(sigmoid_2nd_deriv)
                    infl_dex = sigmoid_2nd_deriv.index(infl_deriv)
                    inflect = sigdatdom[infl_dex]
                    #inflect = poptsi[4]*60.0
                    print 'inflection!', inflect, self.well_id
                    inflectdelts = [abs(t_ - inflect) for t_ in t]
                    inflectdex = inflectdelts.index(min(inflectdelts))
                    try: mi = vals[inflectdex]
                    except IndexError:
                        print 'could not find OD threshold from sigmoidal-fit!',self.well_id
                        mi = self.mi_thresh
                lodex,midex,hidex = self.find_threshold_indices(
                    lo,mi,hi,self.cufit_data_sig.scalars)
                lodex += siglodex
                midex += siglodex
                hidex += siglodex
            except RuntimeError:
                print 'runtimeerror: scipy could not fit a curve for well',self.well_id
                self.cufit_data_sig = ldc.scalars(label = 'sigmoid-fit', 
                    scalars = [], domain = [], override_domain = True)
                lodex,midex,hidex = self.find_threshold_indices(lo,mi,hi,vals)
            except TypeError:
                print 'typeerror: scipy could not fit a curve for well', self.well_id
                self.cufit_data_sig = ldc.scalars(label = 'sigmoid-fit', 
                    scalars = [], domain = [], override_domain = True)
                lodex,midex,hidex = self.find_threshold_indices(lo,mi,hi,vals)

            self.threshlines = []
            self.threshlines.append(ldc.scalars(label = '__skip__lo', 
                scalars = [lo,lo], override_domain = True, color = 'black',  
                linestyle = '--',linewidth = 2.0,
                domain = [t[0],t[len(vals) - 1]]))
            self.threshlines.append(ldc.scalars(label = '__skip__mi', 
                scalars = [mi,mi], override_domain = True, color = 'black',  
                linestyle = '--',linewidth = 2.0,
                domain = [t[0],t[len(vals) - 1]]))
            self.threshlines.append(ldc.scalars(label = '__skip__hi', 
                scalars = [hi,hi], override_domain = True, color = 'black', 
                linestyle = '--',linewidth = 2.0,
                domain = [t[0],t[len(vals) - 1]]))

            expxrelev = t[lodex:midex]/60.0
            expyrelev = vals[lodex:midex]
            if len(expyrelev) == 0:
                print 'empty data set for well', self.well_id
                return 0.0, 0.0
            try:
                poptex, pcovex = cufit(expo, expxrelev, expyrelev)
                self.cufit_data_exp = ldc.scalars(label = 'expo-fit', 
                    scalars = expo(expxrelev, *poptex), domain = 60.0*expxrelev, 
                    override_domain = True, color = 'g', linewidth = 1.5)
            except RuntimeError:
                print 'runtimeerror: scipy could not fit a curve for well', self.well_id
                return 0.0, 0.0
            except TypeError:
                print 'runtimeerror: scipy could not fit a curve for well', self.well_id
                return 0.0, 0.0

            self.lo_thresh = lo
            self.mi_thresh = mi
            self.hi_thresh = hi
            sigfigs = self.significant_figures
            grate = round(poptex[1]/60.0,sigfigs)
            if grate == 0.0:
                print 'growth rate rounded to 0.0!', self.well_id
                return 0.0,0.0
            dbling = round(np.log(2)/grate,sigfigs)
	    return dbling, grate

class optical_density_block(obs_data_block):
	def __init__(self, *args, **kwargs):                      
		self.impose_default('OD_threshold_low', 0.1)
		self.impose_default('OD_threshold_middle', 0.6)
		self.impose_default('OD_threshold_high', 1.0)
                self.impose_default('selected_row',None)
		self.current_tab_index_tablebook = 0
		#lo = 0.05
		#hi = 0.7
		self.impose_default('_well_select_', None)
                kwargs['is_OD_block'] = True             
		obs_data_block.__init__(self, *args, **kwargs)

        def get_threshold_index(self, vals, th):
	    lodelts = [abs(x-th) for x in vals]
            deldex = lodelts.index(min(lodelts))
            return deldex

        def get_lo_hi_threshold_indexes(self, lo, mi, hi, vals = None):
            indices = []
            #for da in self.data.data:
            #    vals = da.scalars
            if vals is None:
                daters = [d.data for d in self.well_mobjs.values()]
                vcnt = len(daters[0].scalars)
                vals = []
                for tdx in xrange(vcnt):
                    coll = [d.scalars[tdx] for d in daters]
                    me = np.mean(coll)
                    vals.append(me)
            lodex = self.get_threshold_index(vals, lo)
            midex = self.get_threshold_index(vals, mi)
            hidex = self.get_threshold_index(vals, hi)
            #cond_cnt = len(self._cond_key_)
            #indices = indices[cond_cnt:]
            return lodex, midex, hidex

        def impose_global_thresholds(self):
            wobjs = self.well_mobjs
            glow = self.OD_threshold_low
            gmid = self.OD_threshold_middle
            ghigh = self.OD_threshold_high
            for ke in wobjs.keys():
                wobj = wobjs[ke]
                wobj.lo_thresh = glow
                wobj.mi_thresh = gmid
                wobj.hi_thresh = ghigh
            self.recalculate_doubling()

	def recalculate_doubling(self):
            for ke in self.well_mobjs.keys():
                self.recalculate_individual_doubling(ke)
            self.recalculate_replicate_doubling()

        def recalculate_individual_doubling(self, well):
            wobj = self.well_mobjs[well]
            tdata = self.cond_mobjs['Time'].data
            wobj.process(tdata)
	    self.rewidget(True)

        def recalculate_replicate_doubling(self):
            for ke in self.replicate_mobjs.keys():
                robj = self.replicate_mobjs[ke]
                robj.process()
                #wobjs = robj.wells
            self.rewidget(True)

        def change_thresh_callback(self, spinwidg):
            if hasattr(spinwidg._modu_inst_, 'well_id'):
                well = spinwidg._modu_inst_.well_id
            else: return
            self.recalculate_individual_doubling(well)
            dtime_txtbox = self.__dict__['_modu_dtime_txtbox_'+well][0]
            grate_txtbox = self.__dict__['_modu_grate_txtbox_'+well][0]
            newdtime = self.well_mobjs[well].dtime
            newgrate = self.well_mobjs[well].grate
            dtime_txtbox.setText(str(newdtime))
            grate_txtbox.setText(str(newgrate))
            if well == self.selected_row: self.redraw_plot()

        def change_table_selection(self, table, rowlabel):
            self.selected_row = rowlabel
            self.redraw_plot(rowlabel)

        def redraw_plot(self, well = None):
            if not well in self._well_key_:
                if self.selected_row is None: return
                else: well = self.selected_row
            qplot = self.qplot[0]
            data = self.get_well_data(well)
            ptype = 'lines'
            qplot.plot(data,'time','OD',
                'OD of ' + well, ptype = ptype)
            print 'want to show plot for well', well
            
        def get_well_data(self, well):
            datacont = lfu.data_container()
            xrelev = self.cond_mobjs['Time'].data
            yrelev = self.well_mobjs[well].data
            cudataexp = self.well_mobjs[well].cufit_data_exp
            cudatasig = self.well_mobjs[well].cufit_data_sig
            threshlines = self.well_mobjs[well].threshlines
            datacont.data = [xrelev,yrelev,cudataexp,cudatasig]+threshlines
            datacont.x_log = False
            datacont.y_log = False
            datacont.active_targs = [well,
                'expo-fit','sigmoid-fit',
                '__skip__lo','__skip__mi','__skip__hi']
            datacont.xdomain = 'Time'
            return datacont

        def calculate_rep_table(self):
            repheads = [
                'Mean Doubling Time', 
                'Mean Growth Rate', 
                'Stddev Of Doubling Time', 
                'Stddev Of Growth Rate']
            reptemps = []
            reprows = self._replicate_key_[:]
            for rdx,ro in enumerate(reprows):
                replic = self.replicate_mobjs[ro]
                #reptemps.append([])
                mdtime = str(replic.mean_doubling_time)
                mgrate = str(replic.mean_growth_rate)
                stddtime = str(replic.stddev_doubling_time)
                stdgrate = str(replic.stddev_growth_rate)
                rotemps = [
                    lgm.interface_template_gui(
                        widgets = ['text'], 
                        read_only = [True], 
                        initials = [[mdtime]]), 
                    lgm.interface_template_gui(
                        widgets = ['text'], 
                        read_only = [True], 
                        initials = [[mgrate]]), 
                    lgm.interface_template_gui(
                        widgets = ['text'], 
                        read_only = [True], 
                        initials = [[stddtime]]), 
                    lgm.interface_template_gui(
                        widgets = ['text'], 
                        read_only = [True], 
                        initials = [[stdgrate]]), 
                    ]
                reptemps.append(rotemps)
            return repheads, reprows, reptemps

	def set_settables(self, *args, **kwargs):
		window = args[0]
		self.handle_widget_inheritance(*args, **kwargs)
		self.widg_templates.append(
                    lgm.interface_template_gui(
		        widgets = ['button_set'], 
			bindings = [[
                            [lgb.create_reset_widgets_wrapper(
                                window, self.recalculate_doubling),
                                    self.redraw_plot], 
                            [lgb.create_reset_widgets_wrapper(
                                window, self.impose_global_thresholds), 
                                    self.redraw_plot]]], 
			labels = [['Recalculate Doubling Times', 
                            'Impose Global Thresholds']]))
                wobjs = self.well_mobjs
                heads = ['Low OD Cutoff', 'Middle OD Cutoff', 
                    'High OD Cutoff', 'Doubling Time', 'Growth Rate']
                rows = ['Global'] + [we for we in self._well_key_]
                temps = []
                for row in rows:
                    rowtemps = []
                    for head in heads:
                        if head == 'Low OD Cutoff':
                            if row == 'Global':
                                inst = self
                                key = 'OD_threshold_low'
                            else:
                                inst = wobjs[row]
                                key = 'lo_thresh'
                            rowtemps.append(
                                lgm.interface_template_gui(
                                    doubles = [[True]], 
                                    widgets = ['spin'], 
                                    single_steps = [[0.01]], 
                                    callbacks = [[self.change_thresh_callback]], 
                                    instances = [[inst]], 
                                    keys = [[key]], 
                                    initials = [[inst.__dict__[key]]], 
                                    ))
                        elif head == 'Middle OD Cutoff':
                            if row == 'Global':
                                inst = self
                                key = 'OD_threshold_middle'
                            else:
                                inst = wobjs[row]
                                key = 'mi_thresh'
                            rowtemps.append(
                                lgm.interface_template_gui(
                                    doubles = [[True]], 
                                    widgets = ['spin'], 
                                    single_steps = [[0.01]], 
                                    callbacks = [[self.change_thresh_callback]], 
                                    instances = [[inst]], 
                                    keys = [[key]], 
                                    initials = [[inst.__dict__[key]]], 
                                    ))
                        elif head == 'High OD Cutoff':
                            if row == 'Global':
                                inst = self
                                key = 'OD_threshold_high'
                            else:
                                inst = wobjs[row]
                                key = 'hi_thresh'
                            rowtemps.append(
                                lgm.interface_template_gui(
                                    doubles = [[True]], 
                                    widgets = ['spin'], 
                                    single_steps = [[0.01]], 
                                    instances = [[inst]], 
                                    keys = [[key]], 
                                    callbacks = [[self.change_thresh_callback]], 
                                    initials = [[inst.__dict__[key]]], 
                                    ))
                        elif head == 'Doubling Time':
                            if not row == 'Global':
                                wobj = wobjs[row]
                                dtime = str(wobj.dtime)
                                rowtemps.append(
                                    lgm.interface_template_gui(
                                        widgets = ['text'], 
                                        read_only = [True], 
                                        initials = [[dtime]], 
                                        handles = [(self, 
                                            '_modu_dtime_txtbox_' + row)], 
                                        ))
                            else: rowtemps.append(None)
                        elif head == 'Growth Rate':
                            if not row == 'Global':
                                wobj = wobjs[row]
                                grate = str(wobj.grate)
                                rowtemps.append(
                                    lgm.interface_template_gui(
                                        widgets = ['text'], 
                                        read_only = [True], 
                                        initials = [[grate]], 
                                        handles = [(self, 
                                            '_modu_grate_txtbox_' + row)], 
                                        ))
                            else: rowtemps.append(None)
                    temps.append(rowtemps)
                table_for_all_wells_template = \
                        lgm.interface_template_gui(
                            layout = 'vertical', 
                            widgets = ['check_set'], 
                            append_instead = [False], 
                            labels = [['Autofind Thresholds Using Sigmoid-Fit']], 
                            instances = [[self]], 
			    #callbacks = [[lgb.create_reset_widgets_wrapper(
			    #    window, self.recalculate_doubling)]], 
                            keys = [['override_thresholds']])
                table_plot_template =\
                        lgm.interface_template_gui(
                                layout = 'horizontal', 
                                widgets = ['table'], 
                                labels = [[heads, rows]],  
                                callbacks = [[self.change_table_selection]], 
                                templates = [temps])
                data = self.get_well_data('A1')
                table_plot_template +=\
                    lgm.interface_template_gui(
                        widgets = ['plot'], 
                        handles = [(self, 'qplot')], 
                        datas = [data])
                table_for_all_wells_template +=\
                    lgm.interface_template_gui(
                        widgets = ['panel'], 
                        templates = [[table_plot_template]])

                repheads, reprows, reptemps = self.calculate_rep_table()
                table_for_reps_template =\
                    lgm.interface_template_gui(
                        layout = 'vertical', 
                        widgets = ['table'], 
                        labels = [[repheads, reprows]], 
                        templates = [reptemps])
                table_for_reps_template +=\
                    lgm.interface_template_gui(
                        widgets = ['button_set'], 
                        bindings = [[lgb.create_reset_widgets_wrapper(
                            window, self.recalculate_replicate_doubling)]], 
                        labels = [['Recalculate Replicate Data']])
                        
                tablepages = [('Raw',[table_for_all_wells_template]), 
                            ('Replicates', [table_for_reps_template])]
		self.widg_templates.append(
                    lgm.interface_template_gui(
			widgets = ['tab_book'], 
			pages = [tablepages], 
			initials = [[self.current_tab_index_tablebook]], 
			instances = [[self]], 
			keys = [['current_tab_index_tablebook']]))

                '''
                self.widg_templates.append(
                        lgm.interface_template_gui(
                            widgets = ['check_set'], 
                            append_instead = [False], 
                            labels = [['Autofind Thresholds Using Sigmoid-Fit']], 
                            instances = [[self]], 
			    #callbacks = [[lgb.create_reset_widgets_wrapper(
			    #    window, self.recalculate_doubling)]], 
                            keys = [['override_thresholds']]))
                self.widg_templates.append(
                        lgm.interface_template_gui(
                                layout = 'horizontal', 
                                widgets = ['table'], 
                                labels = [[heads, rows]],  
                                callbacks = [[self.change_table_selection]], 
                                templates = [temps]))
                data = self.get_well_data('A1')
                self.widg_templates[-1] +=\
                    lgm.interface_template_gui(
                        widgets = ['plot'], 
                        handles = [(self, 'qplot')], 
                        datas = [data])
                '''
		obs_data_block.set_settables(self, *args, from_sub = True)

class data_pool(data_block):
	def __init__(self, *args, **kwargs):
		self.data = []
		data_block.__init__(self, [], **kwargs)
	def append(self, value):
		value.parent = self
		self.data.append(value)
                if value.is_OD_block: self.OD_block_index = len(self.data) - 1
		self._children_ = self.data
	def _read_(self, *args, **kwargs):
		for bl in self.data: bl._read_(*args, **kwargs)
	def deep_parse(self, *args, **kwargs):
		for bl in self.data: bl.deep_parse(*args, **kwargs)
	def get_targets(self):
		return [bl.label for bl in self.data]
        def get_OD_block(self):
                return self.data[self.OD_block_index]

class plate_reader_analyzer(lfu.modular_object_qt):

	def __init__(self, *args, **kwargs):
		self.impose_default('parsed_data',lfu.data_container(),**kwargs)
		self.settings_manager = lset.settings_manager(parent = self, 
		    	filename = 'plate_reader_analyzer_settings.txt')
		self.settings = self.settings_manager.read_settings()
		in_dat = lset.get_setting('default_input_data')
		in_tmp = lset.get_setting('default_input_template')
		self.impose_default('input_data_file',in_dat,**kwargs)
		#self.impose_default('input_tmpl_file',in_tmp,**kwargs)
                self.impose_default('significant_figures',3,**kwargs)
                self.impose_default('reps_per_rep', '3', **kwargs)
                self.impose_default('first_blank_well', None, **kwargs)
		self.current_tab_index = 0
		self.current_tab_index_outputs = 0
		self.postprocess_plan = lpp.post_process_plan(
			label = 'Post Process Plan', parent = self)
		self.postprocess_plan._display_for_children_ = True
		lfu.modular_object_qt.__init__(self, *args, **kwargs)
		self._children_ = [self.postprocess_plan]

	def parse_inputs(self):
		#self.parse_template()
		self.parse_data()
        '''
	def parse_template(self):
		fipath = self.input_tmpl_file
                if not os.path.isfile(fipath):
			print 'cannot find template file:', self.input_tmpl_file
			return
		with open(fipath, 'r') as handle: lines = handle.readlines()
		read = OrderedDict()
		comments = []
		for li in lines:
			l = li.strip()
			if l.startswith('#'): comments.append(l)
			elif not l == '':
				if l.startswith('<') and l.endswith('>'):
					head = l[1:-1].split(',')
					for h in head: read[h] = []
				else:
					body = l.split(',')
					for h, b in zip(read.keys(), body):
						read[h].append(b)
		read['comments'] = comments
		self.template_read = read
        '''

	def parse_data(self):

		def to_blocks(lines):
			blocks = []
			for li in lines:
				if li.startswith('\n') or li.startswith('\r'):
					blocks.append([])
				else:
					try: blocks[-1].append(li[:-1])
					except IndexError: blocks.append([])
			bl_objs = []
			[bl_objs.append(data_block(bl)) for bl in blocks if bl]
			return bl_objs

		fipath = self.input_data_file
		if not os.path.isfile(fipath):
			print 'cannot find data file:', self.input_data_file
			return
		with open(fipath, 'r') as handle: lines = handle.readlines()
		blocks = to_blocks(lines)
		header = header_block(blocks[:3])
		procedure = procedure_block(blocks[3:5])
		layout = layout_block(blocks[5])
		obs_blocks = blocks[6:]
		obs_heads = obs_blocks[::2]
		obs_data = obs_blocks[1::2]
		data = data_pool(parent = self)
		self._children_ = [data]
		_OD = 600
		_OD_lam = 600
		_OD_label = str(_OD_lam)+'nm:'+str(_OD)

                for pair in zip(obs_heads, obs_data):
                        bl_label = swap_check(data_block.swaps, pair[0].raw[0])
			if bl_label == _OD_label:
				data.append(optical_density_block(
						pair, label = bl_label))
			else: data.append(obs_data_block(pair, label = bl_label))

		self.read = {}
		self.read['header'] = header
		self.read['procedure'] = procedure
		self.read['layout'] = layout
		self.read['data'] = data
		for key in self.read.keys(): self.read[key]._read_(self)
		self.read['data'].deep_parse(self)
		daters = self.read['data'].get_targets()
		self.postprocess_plan._always_sourceable_ = daters
		self.postprocess_plan.rewidget(True)
		self.rewidget(True)

	def analyze_data(self):
		try:
			print 'performing analysis...'
			check = time.time()
			self.postprocess_plan(self)
			print 'duration of analysis: ', time.time() - check
			return True
		except:
			traceback.print_exc(file=sys.stdout)
			print 'failed to run post processes'
			return False

	def produce_output(self):
		print 'producing output...'
		check_0 = time.time()
		for dex, dat in enumerate(self.read['data'].data):
			if dat.output.must_output():
				dat.provide_axes_manager_input()
				dat.output(dat.data)
			else: print 'skipping output...', dat.output.label
		print 'produced output: ', time.time() - check_0

	def open_file(self):
		print 'open a file'

	def make_tab_book_pages(self, *args, **kwargs):
		window = args[0]
		inpt_dir = os.getcwd()
		front_page = lgm.interface_template_gui(
				widgets = ['file_name_box'], 
				layout = 'grid', 
                                widg_spans = [(1,2)], 
                                widg_positions = [(2,0)],
				keys = [['input_data_file']], 
				instances = [[self]], 
                                maximum_sizes = [[(50, 200)]], 
				initials = [[self.input_data_file, 
					'Possible Inputs (*.txt)', 
					inpt_dir]], 
				labels = [['Choose Filename']], 
				box_labels = ['Input Data File'])
                '''
		front_page += lgm.interface_template_gui( 
                                widg_positions = [(3,0)],
                                widg_spans = [(1,2)], 
				widgets = ['file_name_box'], 
				keys = [['input_tmpl_file']], 
				instances = [[self]], 
                                maximum_sizes = [[(50, 200)]], 
				initials = [[self.input_tmpl_file, 
			        	'Possible Inputs (*.txt)', 
					inpt_dir]], 
				labels = [['Choose Filename']], 
				box_labels = ['Input Template File'])
                '''
                front_page += lgm.interface_template_gui(
                                widg_positions = [(1,0)],
                                widgets = ['radio'], 
                                labels = [['2', '3', '4']], 
                                box_labels = ['Replicas Per Replicate'], 
                                initials = [[self.reps_per_rep]], 
                                instances = [[self]], 
                                keys = [['reps_per_rep']])
                if hasattr(self, '_well_key_'):
                    wells = self._well_key_
                else: wells = []
                if self.first_blank_well is None and wells:
                    self.first_blank_well = wells[-1*self.reps_per_rep]
                    print 'FIST BLANK WELL', self.first_blank_well
                front_page += lgm.interface_template_gui(
                                widg_positions = [(1,1)],
                                widgets = ['selector'], 
                                instances = [[self]], 
                                keys = [['first_blank_well']], 
                                labels = [wells], 
                                #labels = [self._well_key_], 
                                box_labels = ['First Blank Well'], 
                                initials = [[self.first_blank_well]])
		front_page += lgm.interface_template_gui(
                                widg_positions = [(0,0)], 
				widgets = ['button_set'], 
				layouts = ['horizontal'], 
				bindings = [[lgb.create_reset_widgets_wrapper(
						window, self.parse_inputs), 
					self.analyze_data, self.produce_output]], 
				labels = [['Parse Input File', 
				    'Analyze Parsed Data', 'Produce Output']])
		pages = [('Main', [front_page])]
		output_pages = []
		if hasattr(self, 'read'):
			[output_pages.append((bl.label, 
				bl.output.widg_templates)) 
				for bl in self.read['data'].data]
		for proc in self.postprocess_plan.post_processes:
			try:
				output_pages.append((proc.label, 
					proc.output.widg_templates))
			except AttributeError:
				proc.output.set_settables(*args, **kwargs)
				output_pages.append((proc.label, 
					proc.output.widg_templates))
		output_tabs = lgm.interface_template_gui(
			widgets = ['tab_book'], 
			pages = [output_pages], 
			initials = [[self.current_tab_index_outputs]], 
			instances = [[self]], 
			keys = [['current_tab_index_outputs']])
		pages.append(('Output', [output_tabs]))
		self.postprocess_plan.set_settables(window, self)
		pp_page = lgm.interface_template_gui(
			widgets = ['panel'], 
			templates = [self.postprocess_plan.widg_templates])
		pages.append(('Post Procecssing', [pp_page]))
		if hasattr(self, 'read'):
			[pages.append((bl.label, bl.widg_templates)) 
						for bl in self.read['data'].data]
		return pages

	def rewidget__children_(self, *args, **kwargs):
		kwargs['infos'] = (kwargs['infos'][0], self)
		for child in self._children_:
			if child.rewidget(**kwargs):
				child.set_settables(*kwargs['infos'])

	def set_settables(self, *args, **kwargs):
		window = args[0]
		self.handle_widget_inheritance(*args, **kwargs)
		#gear_icon_path = os.path.join(
		#	os.getcwd(), 'resources', 'gear.png')
		wrench_icon_path = lfu.get_resource_path('wrench_icon.png')
		center_icon_path = lfu.get_resource_path('center.png')
		refresh_icon_path = lfu.get_resource_path('refresh.png')
		open_icon_path = lfu.get_resource_path('open.png')
		quit_icon_path = lfu.get_resource_path('quit.png')
		settings_ = lgb.create_action(parent = window, label = 'Settings', 
			bindings = lgb.create_reset_widgets_wrapper(
			window, self.change_settings), icon = wrench_icon_path, 
			shortcut = 'Ctrl+Shift+S', statustip = 'Settings')
		self.refresh_ = lgb.create_reset_widgets_function(window)
		update_gui_ = lgb.create_action(parent = window, 
			label = 'Refresh GUI', icon = refresh_icon_path, 
			shortcut = 'Ctrl+G', bindings = self.refresh_, 
			statustip = 'Refresh the GUI (Ctrl+G)')
		open_file = lgb.create_action(parent = window, label = 'Open', 
			bindings = lgb.create_reset_widgets_wrapper(
			window, self.open_file), icon = open_icon_path, 
			shortcut = 'Ctrl+O', statustip = 'Open Input File')
		quit_ = lgb.create_action(parent = window, label = 'Quit', 
			icon = quit_icon_path, shortcut = 'Ctrl+Q', 
				statustip = 'Quit the Application', 
					bindings = window.on_close)
		center_ = lgb.create_action(parent = window, label = 'Center', 
			icon = center_icon_path, shortcut = 'Ctrl+Shift+C', 
					statustip = 'Center Window', 
					bindings = [window.on_resize, 
					window.on_center])
		self.menu_templates.append(
			lgm.interface_template_gui(
				menu_labels = ['&File', '&File', 
					'&File', '&File', '&File'], 
				menu_actions = [settings_, center_, 
					update_gui_, open_file, quit_]))
		self.tool_templates.append(
			lgm.interface_template_gui(
				tool_labels = ['&Tools', '&Tools', 
					'&Tools', '&Tools', '&Tools'], 
				tool_actions = [settings_, center_, 
					update_gui_, open_file, quit_]))
		self.widg_templates.append(
			lgm.interface_template_gui(
				widgets = ['tab_book'], 
				verbosities = [0], 
				pages = [self.make_tab_book_pages(*args, **kwargs)], 
				initials = [[self.current_tab_index]], 
				handles = [(self, 'tab_ref')], 
				instances = [[self]], 
				keys = [['current_tab_index']]))
		lfu.modular_object_qt.set_settables(
			self, *args, from_sub = True)




