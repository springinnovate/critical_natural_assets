"""
Script to mask out the carbon by CNA and do other zonal stats, calculations as:

"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLURFDiff.TOTECOSYSC.2050.nc" * "D:\repositories\critical_natural_assets\overlaps_analysis\RF_diff_Hist.tif"
"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLURSDiff.TOTECOSYSC.2050.nc" * "D:\repositories\critical_natural_assets\overlaps_analysis\RFAFRS_diff_RFAF.tif"
"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLUAFDiff.TOTECOSYSC.2050.nc" * "D:\repositories\critical_natural_assets\overlaps_analysis\RFAF_diff_RF.tif"
"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLURFDiff.TOTECOSYSC.2050.nc" * "D:\repositories\critical_natural_assets\overlaps_analysis\masked_RF_diff_Hist.tif"
"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLURSDiff.TOTECOSYSC.2050.nc" * "D:\repositories\critical_natural_assets\overlaps_analysis\masked_RFAFRS_diff_RFAF.tif"
"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLUAFDiff.TOTECOSYSC.2050.nc" * "D:\repositories\critical_natural_assets\overlaps_analysis\masked_RFAF_diff_RF.tif"

"D:\repositories\critical_natural_assets\data\carbon\sq_km_per_cell.tif"

"""
import os

import numpy
from ecoshard import geoprocessing
from osgeo import gdal

CARBON_RASTERS = {
    'NoLURFDiff': r"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLURFDiff.TOTECOSYSC.2050.tif",
    'NoLURSDiff': r"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLURSDiff.TOTECOSYSC.2050.tif",
    'NoLUAFDiff': r"D:\repositories\critical_natural_assets\data\carbon\SSP126.NoLUAFDiff.TOTECOSYSC.2050.tif",
}

COVERAGE_RASTERS = {
    'RF_diff_Hist': r"D:\repositories\critical_natural_assets\overlaps_analysis\RF_diff_Hist.tif",
    'RFAFRS_diff_RFAF': r"D:\repositories\critical_natural_assets\overlaps_analysis\RFAFRS_diff_RFAF.tif",
    'RFAF_diff_RF': r"D:\repositories\critical_natural_assets\overlaps_analysis\RFAF_diff_RF.tif",
    'cna_RF_diff_Hist': r"D:\repositories\critical_natural_assets\overlaps_analysis\masked_RF_diff_Hist.tif",
    'cna_RFAFRS_diff_RFAF': r"D:\repositories\critical_natural_assets\overlaps_analysis\masked_RFAFRS_diff_RFAF.tif",
    'cna_RFAF_diff_RF': r"D:\repositories\critical_natural_assets\overlaps_analysis\masked_RFAF_diff_RF.tif",
    }

CELL_AREA_KM2_PATH = r"D:\repositories\critical_natural_assets\data\carbon\sq_km_per_cell.tif"

OPERATIONS = [
    ('NoLURFDiff', 'RF_diff_Hist', 'carbon_RF'),
    ('NoLURFDiff', 'cna_RF_diff_Hist', 'carbon_RF_CNA'),
    ('NoLUAFDiff', 'RFAF_diff_RF', 'carbon_RFAF'),
    ('NoLUAFDiff', 'cna_RFAF_diff_RF', 'carbon_RFAF_CNA'),
    ('NoLURSDiff', 'RFAFRS_diff_RFAF', 'carbon_RFAFRS'),
    ('NoLURSDiff', 'cna_RFAFRS_diff_RFAF', 'carbon_RFAFRS_CNA'),
    ]


WORKSPACE_DIR = 'are_cna_carbon_analysis_workspace'
os.makedirs(WORKSPACE_DIR, exist_ok=True)

NODATA_THRESHOLD = -1e30  # anything less than this is nodata lets say


def _calc_carbon_per_cell(
        carbon_density_g_per_sq_m, coverage_pct, sq_km_per_cell,
        carbon_density_nodata, coverage_nodata, target_nodata,
        target_cell_size):
    result = (
        carbon_density_g_per_sq_m * 1000 *  # convert to kg/km^2
        coverage_pct / 100 * sq_km_per_cell *
        # convert area from 0.25 degree cell to smaller
        (target_cell_size/0.25)**2)
    nodata_mask = (
        numpy.isclose(carbon_density_g_per_sq_m, carbon_density_nodata) |
        numpy.isclose(coverage_pct, coverage_nodata))

    result[nodata_mask] = target_nodata
    return result


def main():
    for carbon_density_key, coverage_key, target_prefix in OPERATIONS:
        print(f'processing {target_prefix}')
        target_path = os.path.join(WORKSPACE_DIR, f'{target_prefix}.tif')

        raster_info = geoprocessing.get_raster_info(COVERAGE_RASTERS[coverage_key])

        intermediate_dir = os.path.join(WORKSPACE_DIR, target_prefix)
        os.makedirs(intermediate_dir, exist_ok=True)
        base_raster_path_list = [
            CARBON_RASTERS[carbon_density_key],
            COVERAGE_RASTERS[coverage_key],
            CELL_AREA_KM2_PATH,
        ]

        bounding_box = geoprocessing.merge_bounding_box_list(
            [geoprocessing.get_raster_info(path)['bounding_box']
             for path in base_raster_path_list], 'union')
        aligned_raster_list = [
            os.path.join(intermediate_dir, os.path.basename(path))
            for path in base_raster_path_list
        ]
        for base_raster_path, aligned_raster_path in zip(
                base_raster_path_list, aligned_raster_list):
            geoprocessing.warp_raster(
                base_raster_path, raster_info['pixel_size'],
                aligned_raster_path, 'near', target_bb=bounding_box,
                target_projection_wkt=raster_info['projection_wkt'])
        geoprocessing.align_and_resize_raster_stack(
            base_raster_path_list, aligned_raster_list, ['near']*3,
            raster_info['pixel_size'], 'intersection')
        carbon_density_nodata = geoprocessing.get_raster_info(
            CARBON_RASTERS[carbon_density_key])['nodata'][0]
        coverage_nodata = geoprocessing.get_raster_info(
            COVERAGE_RASTERS[coverage_key])['nodata'][0]
        geoprocessing.raster_calculator(
            [(path, 1) for path in aligned_raster_list] +
            [(carbon_density_nodata, 'raw'), (coverage_nodata, 'raw'),
             (NODATA_THRESHOLD, 'raw'),
             (raster_info['pixel_size'][0], 'raw')],
            _calc_carbon_per_cell,
            target_path,
            gdal.GDT_Float32,
            NODATA_THRESHOLD,
            allow_different_blocksize=True)
        print(f'done with {target_path}')
    print('all done!')


if __name__ == '__main__':
    main()
