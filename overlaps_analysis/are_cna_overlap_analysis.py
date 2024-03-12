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
OVERLAP_THRESHOLD = 3  # Represents the "90%" target for CNA

COUNTRY_VECTOR_PATH = "../data/countries_iso3_md5_6fb2431e911401992e6e56ddf0a9bcda.gpkg"

WORKING_DIR = 'are_cna_overlap_analysis_workspace'
os.makedirs(WORKING_DIR, exist_ok=True)

def main():
    table_file = open('summary.csv', 'w')
    table_file.write(
        'country name,'
        'RF_country_coverage,AF_country_coverage,RS_country_coverage,'
        'RF_country_avg_coverage,AF_country_avg_coverage,RS_country_avg_coverage,'
        'RF_country_cna_coverage,AF_country_cna_coverage,RS_country_cna_coverage,'
        'RF_country_cna_avg_coverage,AF_country_cna_avg_coverage,RS_country_cna_avg_coverage,'
        'CNA_country_percent_cover,'
        'area_sq_km,'
        '\n')
    country_stat_results = collections.defaultdict(lambda: collections.defaultdict(dict))
    task_graph = taskgraph.TaskGraph(WORKING_DIR, -1)
    country_stats_list = []
    for source_index, diff_index in [(1, 0), (2, 1), (3, 2)]:
        source_var = VARIABLES[source_index]
        diff_var = VARIABLES[diff_index]
        source_path = FOREST_COVER_RASTER_MAP[source_var]
        diff_path = FOREST_COVER_RASTER_MAP[diff_var]
        new_forest_raster_path = f'{source_var}_diff_{diff_var}.tif'
        print(f'processing {new_forest_raster_path}')
        diff_task = task_graph.add_task(
            func=geoprocessing.raster_calculator,
            args=(
                [(source_path, 1), (diff_path, 1)], lambda x, y: x-y, new_forest_raster_path,
                gdal.GDT_Float32, 0),
            target_path_list=[new_forest_raster_path],
            kwargs={'allow_different_blocksize': True})
        overlap_raster_info = geoprocessing.get_raster_info(OVERLAP_RASTER_PATH)
        projected_raster_path = os.path.join(WORKING_DIR, f'diff_{os.path.basename(new_forest_raster_path)}')
        warp_task = task_graph.add_task(
            func=geoprocessing.warp_raster,
            args=(
                new_forest_raster_path, overlap_raster_info['pixel_size'], projected_raster_path,
                'near'),
            kwargs={
                'target_bb': overlap_raster_info['bounding_box'],
                'target_projection_wkt': overlap_raster_info['projection_wkt']},
            dependent_task_list=[diff_task],
            target_path_list=[projected_raster_path])

        masked_forest_raster_path = f'masked_{os.path.basename(new_forest_raster_path)}'
        mask_task = task_graph.add_task(
            func=geoprocessing.raster_calculator,
            args=(
                [(projected_raster_path, 1), (OVERLAP_RASTER_PATH, 1)], _mask_raster_op, masked_forest_raster_path,
                gdal.GDT_Float32, 0),
            kwargs={'allow_different_blocksize': True},
            target_path_list=[masked_forest_raster_path],
            dependent_task_list=[warp_task])

        cna_country_stats = task_graph.add_task(
            func=geoprocessing.zonal_statistics,
            args=(
                (masked_forest_raster_path, 1), COUNTRY_VECTOR_PATH),
            kwargs={'working_dir': WORKING_DIR},
            dependent_task_list=[mask_task],
            store_result=True)

        raw_tree_country_stats = task_graph.add_task(
            func=geoprocessing.zonal_statistics,
            args=(
                (projected_raster_path, 1), COUNTRY_VECTOR_PATH),
            kwargs={'working_dir': WORKING_DIR},
            dependent_task_list=[mask_task],
            store_result=True)
        country_stats_list.append((cna_country_stats, raw_tree_country_stats, source_var))

    cna_mask_raster_path = os.path.join(WORKING_DIR, 'cna_mask.tif')
    cna_mask_task = task_graph.add_task(
        func=geoprocessing.raster_calculator,
        args=(
            [(OVERLAP_RASTER_PATH, 1)], lambda x: x>2, cna_mask_raster_path, gdal.GDT_Float32, 0),
        kwargs={'allow_different_blocksize': True},
        target_path_list=[cna_mask_raster_path])
    cna_in_country_stats = task_graph.add_task(
        func=geoprocessing.zonal_statistics,
        args=(
            (cna_mask_raster_path, 1), COUNTRY_VECTOR_PATH),
        kwargs={'working_dir': WORKING_DIR},
        dependent_task_list=[cna_mask_task],
        store_result=True)

    global_stats = collections.defaultdict(lambda: collections.defaultdict(float))
    country_stat_results['0_GLOBAL']['area_sq_km'] = 0
    for cna_country_stats, raw_tree_country_stats, source_var in country_stats_list:
        country_vector = gdal.OpenEx(COUNTRY_VECTOR_PATH, gdal.OF_VECTOR)
        country_layer = country_vector.GetLayer()
        for country_feature in country_layer:
            country_name = country_feature.GetField('nev_name')
            cna_local_stats = cna_country_stats.get()[country_feature.GetFID()]
            global_stats[source_var]['sum'] += cna_local_stats['sum']
            global_stats[source_var]['count'] += cna_local_stats['count']
            global_stats[source_var]['nodata_count'] += cna_local_stats['nodata_count']

            try:
                country_stat_results[country_name][source_var]['cna_avg'] = cna_local_stats['sum']/cna_local_stats['count']
                country_stat_results[country_name][source_var]['cna_coverage'] = 100*cna_local_stats['count']/(cna_local_stats['nodata_count']+cna_local_stats['count'])
            except ZeroDivisionError:
                country_stat_results[country_name][source_var]['cna_avg'] = 0
                country_stat_results[country_name][source_var]['cna_coverage'] = 0


            raw_tree_local_stats = raw_tree_country_stats.get()[country_feature.GetFID()]
            global_stats[source_var]['raw_sum'] += raw_tree_local_stats['sum']
            global_stats[source_var]['raw_count'] += raw_tree_local_stats['count']
            global_stats[source_var]['raw_nodata_count'] += raw_tree_local_stats['nodata_count']
            try:
                country_stat_results[country_name][source_var]['raw_tree_avg'] = raw_tree_local_stats['sum']/raw_tree_local_stats['count']
                country_stat_results[country_name][source_var]['raw_tree_coverage'] = 100*raw_tree_local_stats['count']/(raw_tree_local_stats['nodata_count']+raw_tree_local_stats['count'])
            except ZeroDivisionError:
                country_stat_results[country_name][source_var]['raw_tree_avg'] = 0
                country_stat_results[country_name][source_var]['raw_tree_coverage'] = 0

    country_layer.ResetReading()
    global_cna_count = 0
    global_cna_nodata_count = 0
    for country_feature in country_layer:
        country_name = country_feature.GetField('nev_name')
        country_stat_results[country_name]['area_sq_km'] = country_feature.GetField('area')/1000**2
        cna_stats = cna_in_country_stats.get()[country_feature.GetFID()]
        country_stat_results[country_name]['CNA_country_percent_cover'] = 100*cna_stats['count']/(cna_stats['count']+cna_stats['nodata_count'])
        global_cna_count += cna_stats['count']
        global_cna_nodata_count += cna_stats['nodata_count']
        country_stat_results['0_GLOBAL']['area_sq_km'] += country_stat_results[country_name]['area_sq_km']

    for source_var in global_stats:
        country_stat_results['0_GLOBAL'][source_var]['cna_avg'] = global_stats[source_var]['sum']/global_stats[source_var]['count']
        country_stat_results['0_GLOBAL'][source_var]['cna_coverage'] = global_stats[source_var]['sum']/(global_stats[source_var]['nodata_count']+global_stats[source_var]['count'])
        country_stat_results['0_GLOBAL'][source_var]['raw_tree_avg'] = global_stats[source_var]['raw_sum']/global_stats[source_var]['raw_count']
        country_stat_results['0_GLOBAL'][source_var]['raw_tree_coverage'] = global_stats[source_var]['raw_sum']/(global_stats[source_var]['raw_nodata_count']+global_stats[source_var]['count'])

    country_stat_results['0_GLOBAL']['CNA_country_percent_cover'] = 100*global_cna_count/(global_cna_count+global_cna_nodata_count)
    # what fraction of the pixels have a positive change in diff
    # what's the average of the % of forest pixels in the diff
    # what percent of the country that is CNA overlaps with a positive diff pixel
    # what's the average of the % of forest pixels CNA

    for country_name in sorted(country_stat_results):
        local_stats = country_stat_results[country_name]
        table_file.write(
            f"{country_name},"
            f"{local_stats['RF']['raw_tree_coverage']},{local_stats['RFAF']['raw_tree_coverage']},{local_stats['RFAFRS']['raw_tree_coverage']},"
            f"{local_stats['RF']['raw_tree_avg']},{local_stats['RFAF']['raw_tree_avg']},{local_stats['RFAFRS']['raw_tree_avg']},"
            f"{local_stats['RF']['cna_coverage']},{local_stats['RFAF']['cna_coverage']},{local_stats['RFAFRS']['cna_coverage']},"
            f"{local_stats['RF']['cna_avg']},{local_stats['RFAF']['cna_avg']},{local_stats['RFAFRS']['cna_avg']},"
            f"{country_stat_results[country_name]['CNA_country_percent_cover']},"
            f"{local_stats['area_sq_km']}"
            "\n")

        # 'country name,'
        # 'RF_country_coverage,RFAF_country_coverage,RFAFRS_country_coverage,'
        # 'RF_country_avg_coverage,RFAF_country_avg_coverage,RFAFRS_country_avg_coverage,'
        # 'RF_country_cna_coverage,RFAF_country_cna_coverage,RFAFRS_country_cna_coverage,'
        # 'RF_country_cna_avg_coverage,RFAF_country_cna_avg_coverage,RFAFRS_country_cna_avg_coverage,'
        # '\n')

def _mask_raster_op(base_array, mask_array):
    result = base_array.copy()
    result[mask_array <= OVERLAP_THRESHOLD] = 0
    return result

if __name__ == '__main__':
    main()
