<?php
$gene_file = $_POST['gene_file'];
$cutoff = $_POST['cutoff'];
$top_cutoff = $_POST['top_cutoff'];
$go_file = $_POST['go_file'];
$kind_of_list = $_POST['kind_of_list'];
$output = exec("python hist_cc_web.py $gene_file $cutoff $top_cutoff $go_file $kind_of_list");
echo $output;
?>