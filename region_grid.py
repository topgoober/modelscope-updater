"""Bounded regional grids shared by the job API and isolated decoder."""
import json, math, os

def regional_grid(latitude, longitude):
    latitude, longitude = float(latitude), float(longitude)
    if not math.isfinite(latitude) or not math.isfinite(longitude) or not (24 <= latitude <= 50 and -125 <= longitude <= -65):
        raise ValueError('Select a location within the contiguous-US map')
    # Quantization lets nearby clicks reuse a bounded cache; ~3.3 km N-S sampling.
    latitude, longitude = round(latitude * 2) / 2, round(longitude * 2) / 2
    west = max(-125, min(-69.8, longitude - 2.4))
    north = min(50, max(27.6, latitude + 1.8))
    return dict(west=round(west, 2), north=round(north, 2), step=.03, width=161, height=121)

def display_grids():
    center = os.environ.get('MODELSCOPE_REGION_CENTER')
    if center:
        lat, lon = json.loads(center)
        grid = regional_grid(lat, lon)
        return grid, dict(grid), 1
    return (dict(west=-125., north=50., step=.4, width=149, height=66),
            dict(west=-125., north=50., step=.8, width=75, height=33), 2)
