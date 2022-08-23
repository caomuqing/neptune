import os
import rospy
import rospkg

from qt_gui.plugin import Plugin
from python_qt_binding import loadUi
from python_qt_binding.QtWidgets import QWidget

from behavior_selector.srv import MissionModeChange

# Must match values in MissionModeChange
START = 1
END   = 2
KILL  = 3

class MissionModePlugin(Plugin):

    def __init__(self, context):
        super(MissionModePlugin, self).__init__(context)
        # Give QObjects reasonable names
        self.setObjectName('MissionModePlugin')  

        # Create QWidget
        self._widget = QWidget()
        rp = rospkg.RosPack()
        # Get path to UI file which should be in the "resource" folder of this package
        ui_file = os.path.join(rp.get_path('behavior_selector'), 'resource', 'MissionModePlugin.ui')
        # Extend the widget with all attributes and children from UI file
        loadUi(ui_file, self._widget)
        # Give QObjects reasonable names
        self._widget.setObjectName('MissionModePluginUi')
        # Show _widget.windowTitle on left-top of each plugin (when 
        # it's set in _widget). This is useful when you open multiple 
        # plugins at once. Also if you open multiple instances of your 
        # plugin at once, these lines add number to make it easy to 
        # tell from pane to pane.
        if context.serial_number() > 1:
            self._widget.setWindowTitle(self._widget.windowTitle() + (' (%d)' % context.serial_number()))
        # Add widget to the user interface
        context.add_widget(self._widget)

        self._widget.start_push_button.pressed.connect(self._on_start_pressed)
        self._widget.end_push_button.pressed.connect(self._on_end_pressed)
        self._widget.stop_push_button.pressed.connect(self._on_stop_pressed)

    def _on_start_pressed(self):
        mode = START
        self._change_mode(mode)

    def _on_end_pressed(self):
        mode = END
        self._change_mode(mode)

    def _on_stop_pressed(self):
        mode = KILL
        self._change_mode(mode)

    def _change_mode(self,mode):
        rospy.wait_for_service('change_mode',1)
        try:
            mm = rospy.ServiceProxy('change_mode', MissionModeChange)
            mm(mode)
        except rospy.ServiceException, e:
            print "Service call failed: %s"%e

    
