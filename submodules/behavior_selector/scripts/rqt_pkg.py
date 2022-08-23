#!/usr/bin/env python
import sys

from behaviour_selector.button_module import MissionModePlugin
from rqt_gui.main import Main

plugin = 'rqt_pkg'
main = Main(filename=plugin)
sys.exit(main.main(standalone=plugin))
