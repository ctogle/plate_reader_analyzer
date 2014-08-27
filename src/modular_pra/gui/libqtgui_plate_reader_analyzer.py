import modular_core.libfundamental as lfu
import modular_core.gui.libqtgui as lqg
lgm = lqg.lgm
lgd = lqg.lgd
lgb = lqg.lgb
import modular_pra.libplate_reader_analyzer as lprog

import os

class application_plate_reader_analyzer(lqg.application):
	_content_ = [lprog.plate_reader_analyzer()]
	gear_icon = lfu.get_resource_path('gear.png')

	def __init__(self, *args, **kwargs):
		lqg.application.__init__(self, *args, **kwargs)
		x, y = lfu.convert_pixel_space(1024, 256)
		x_size, y_size = lfu.convert_pixel_space(512, 512)
		self._standards_ = {
			'title' : 'plate_reader_analyzer', 
			'geometry' : (x, y, x_size, y_size), 
			'window_icon' : self.gear_icon}
		lqg._window_.apply_standards(self._standards_)
		lqg.application.setStyle(lgb.create_style('plastique'))

_application_ = application_plate_reader_analyzer
_application_locked_ = False

