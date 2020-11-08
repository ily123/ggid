"""This module generates a color gradient."""
import matplotlib.cm
import matplotlib.colors
import numpy as np

class ColorGradientGenerator:
    """Generates color gradient."""

    def __init__(self, lower_bound = 0, upper_bound = 10):
        """Inits class with lower and upper bound for gradient."""

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def create_color_map(self, palette = 'Blues'):
        """Sets colormap."""

        num_bins = 10
        color_map = matplotlib.cm.get_cmap(palette, num_bins)
        self.scalar_map = matplotlib.cm.ScalarMappable(cmap=color_map)

    def create_color_map2(self, base_color = None):
        if not base_color:
            base_color = [255, 0, 0] # red
        N = 256
        vals = np.ones((N, 4))
        vals[:, 0] = np.linspace(1, base_color[0]/256, N)
        vals[:, 1] = np.linspace(1, base_color[1]/256, N)
        vals[:, 2] = np.linspace(1, base_color[2]/256, N)
        color_map = matplotlib.colors.ListedColormap(vals)
        self.scalar_map = matplotlib.cm.ScalarMappable(cmap=color_map)

    def map_colors(self, scores):
        """Maps scores to colors in gradient."""

        rgba_colors = self.scalar_map.to_rgba(scores)
        hex_colors = [matplotlib.colors.to_hex(clr) for clr in rgba_colors]
        return hex_colors
