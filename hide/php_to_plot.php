<?php
$extension=".txt";
switch ($_POST["data_choice"])
{
case 1:
  $dataset = "bkz_gene_list";
  $extension=".csv";
  break;
case 2:
  $dataset = "3-cell-hep-d12_versus_2-cell-hep-d12";
  break;
case 3:
  $dataset = "3-cell-hep-d12_versus_CS-d12";
  break;
case 4:
  $dataset = "3-cell-hep-d12_versus_HM-d12";
  break;
case 5:
  $dataset = "2-cell-hep-d12_versus_CS-d12";
  break;
case 6:
  $dataset = "2-cell-hep-d12_versus_HM-d12";
  break;
case 7:
  $dataset = "CS-d12_versus_HM-d12";
  break;
default:
  $dataset = "3-cell-hep-d12_versus_2-cell-hep-d12";
}
switch ($_POST["onto_choice"])
{
case 1:
  $go_file = "human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv";
  break;
case 2:
  $go_file = "msigdb.txt";
  break;
default:
  $go_file = "human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv";
}
$dataset = $dataset.$extension;

$cutoff = $_POST['cutoff'];
$top_cutoff = $_POST['top_cutoff'];
$kind_of_list = $_POST['kind_of_list'];
$file = $_POST['python'];
$output = exec("python $file $dataset $cutoff $top_cutoff $go_file $kind_of_list");
echo $output;
?>