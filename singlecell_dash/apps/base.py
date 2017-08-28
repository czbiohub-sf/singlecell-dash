"""Building blocks for writing modular callbacks for plotly"""
import dash_html_components as html

CONFIG_DICT = {'modeBarButtonsToRemove': ['sendDataToCloud',
                                          'pan2d',
                                          'zoomIn2d',
                                          'zoomOut2d',
                                          'autoScale2d',
                                          'resetScale2d',
                                          'hoverCompareCartesian',
                                          'hoverClosestCartesian',
                                          'toggleSpikelines'],
               'displaylogo': False
               }


class BaseBlock:

    config_dict = CONFIG_DICT.copy()

    def __init__(self, app=None):
        self.app = app

        if self.app is not None and hasattr(self, 'callbacks'):
            self.callbacks(self.app)

    @property
    def layout(self):
        return html.Div()
