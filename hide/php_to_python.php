<?php
$gene_file = $_POST['gene_file'];
$min_level = $_POST['min_level'];
$cutoff = $_POST['cutoff'];
$top_cutoff = $_POST['top_cutoff'];
$kind = $_POST['k'];
$go_file = $_POST['go_file'];
$kind_of_list = $_POST['kind_of_list'];
$output = exec("python test_ajax5.py $gene_file $min_level $cutoff $top_cutoff $kind $go_file $kind_of_list");
$split_output = explode(",", $output);
echo $split_output[0].','.$split_output[1].','.$split_output[2];
?>