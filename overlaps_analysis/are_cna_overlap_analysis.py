"""
Script that measures the overlap between Stephanie's reforestation
and our old critical natural assets.

Reforestation-Historical
(Reforestation+Aforestation)-Historical
(Reforestation+Aforestation+Restoration)-Historical

CLM5.REFOREST.Hist
CLM5.REFOREST.RF
CLM5.REFOREST.RFAF
CLM5.REFOREST.RFAFRS


"data/CLM5 Deg025 Raw Data-20240130T174114Z-001/CLM5 Deg025 Raw Data/CLM5.REFOREST.RF.001/CLM5.REFOREST.RF.001.2050.PCT_TREE.flt"
"""
import collections
import os

from ecoshard import geoprocessing
from ecoshard import taskgraph
from osgeo import gdal

FOREST_COVER_PATH_PATTERN = "../data/CLM5 Deg025 Raw Data-20240130T174114Z-001/CLM5 Deg025 Raw Data/CLM5.REFOREST.{var}.001/CLM5.REFOREST.{var}.001.2050.PCT_TREE.flt"

FOREST_COVER_RASTER_MAP = {
    'Hist': "../data/CLM5 Deg025 Raw Data-20240130T174114Z-001/CLM5 Deg025 Raw Data/CLM5.REFOREST.Hist.001/CLM5.REFOREST.Hist.001.2015.PCT_TREE.flt",
    'RF': "../data/CLM5 Deg025 Raw Data-20240130T174114Z-001/CLM5 Deg025 Raw Data/CLM5.REFOREST.RF.001/CLM5.REFOREST.RF.001.2050.PCT_TREE.flt",
    'RFAF': "../data/CLM5 Deg025 Raw Data-20240130T174114Z-001/CLM5 Deg025 Raw Data/CLM5.REFOREST.RFAF.001/CLM5.REFOREST.RFAF.001.2050.PCT_TREE.flt",
    'RFAFRS': "../data/CLM5 Deg025 Raw Data-20240130T174114Z-001/CLM5 Deg025 Raw Data/CLM5.REFOREST.RFAFRS.001/CLM5.REFOREST.RFAFRS.001.2050.PCT_TREE.flt"
}

VARIABLES = [
    'Hist',
    'RF',
    'RFAF',
    'RFAFRS']

OVERLAP_RASTER_PATH = "../data/local_NCP_all_targets/local_NCP_land_all_targets_md5_7ccece.tif"
OVERLAP_THRESHOLD = 3

COUNTRY_VECTOR_PATH = "../data/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg"

WORKING_DIR = 'are_cna_overlap_analysis_workspace'
os.makedirs(WORKING_DIR, exist_ok=True)

def main():
    table_file = open('summary.csv', 'w')
    table_file.write('country name,RF,RFAF,RFAFRS\n')
    country_stat_results = collections.default_dict(dict)
    task_graph = taskgraph.TaskGraph(WORKING_DIR, -1)
    for source_index, diff_index in [(1, 0), (2, 0), (3, 0)]:
        source_var = VARIABLES[source_index]
        diff_var = VARIABLES[diff_index]
        source_path = FOREST_COVER_RASTER_MAP[source_var]
        diff_path = FOREST_COVER_RASTER_MAP[diff_var]
        new_forest_raster_path = f'{source_var}_diff_{diff_var}.tif'
        print(f'processing {new_forest_raster_path}')
        geoprocessing.raster_calculator(
            [(source_path, 1), (diff_path, 1)], lambda x, y: x-y, new_forest_raster_path,
            gdal.GDT_Float32, None, allow_different_blocksize=True)
        overlap_raster_info = geoprocessing.get_raster_info(OVERLAP_RASTER_PATH)
        projected_raster_path = os.path.join(WORKING_DIR, f'diff_{os.path.basename(new_forest_raster_path)}')
        geoprocessing.warp_raster(
            new_forest_raster_path, overlap_raster_info['pixel_size'], projected_raster_path,
            'near', target_bb=overlap_raster_info['bounding_box'],
            target_projection_wkt=overlap_raster_info['projection_wkt'])

        masked_forest_raster_path = f'masked_{os.path.basename(new_forest_raster_path)}'
        geoprocessing.raster_calculator(
            [(projected_raster_path, 1), (OVERLAP_RASTER_PATH, 1)], _mask_raster_op, masked_forest_raster_path,
            gdal.GDT_Float32, None, allow_different_blocksize=True)

        country_stats = task_graph.add_task(
            func=geoprocessing.zonal_statistics,
            args=(
                (masked_forest_raster_path, 1), COUNTRY_VECTOR_PATH),
            kwargs={'working_dir': WORKING_DIR},
            store_result=True)
        country_vector = gdal.OpenEx(COUNTRY_VECTOR_PATH, gdal.OF_VECTOR)
        country_layer = country_vector.GetLayer()
        for country_feature in country_layer:
            country_name = country_feature.GetField('nev_name')
            local_stats = country_stats.get()[country_feature.GetFID()]
            avg_pct = local_stats['sum']/local_stats['count']
            country_stat_results[country_name][diff_var] = avg_pct

    for country_name in sorted(country_stats):
        local_stats = country_stats[country_name]
        table_file.write(f"{country_name},{local_stats['RF']},{local_stats['RFAF']},{local_stats['RFAFRS']}\n")

def _mask_raster_op(base_array, mask_array):
    result = base_array.copy()
    result[mask_array <= OVERLAP_THRESHOLD] = 0
    return result

if __name__ == '__main__':
    main()
