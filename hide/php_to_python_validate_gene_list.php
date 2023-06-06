<?php


$help = "";

$extension=".txt";

$allowedExts = array("txt");
$temp = explode(".", $_FILES["file_gene_list"]["name"]);
$file_extension = end($temp);

if (($_FILES["file_gene_list"]["type"] == "text/plain")
   && ($_FILES["file_gene_list"]["size"] < 2000000)
   && in_array($file_extension, $allowedExts)) {
    if ($_FILES["file_gene_list"]["error"] > 0) {
        $help .="<li>Error with the data file: ".$_FILES["file_gene_list"]["error"]."</li>";
    } else {
        $uploadfile = basename($_FILES['file_gene_list']['name']);
        if (move_uploaded_file($_FILES['file_gene_list']['tmp_name'], "uploaded_data/".$uploadfile)) {
            $dataset = "uploaded_data/".$temp[0];
        } 
    }
} else {
    $help .= "<li>Invalid file for data. Check that you chose a txt file.</li>";
}

$dataset = $dataset.$extension; 

switch ($_POST["onto_choice"])
{
case 1:
  $go_file = "human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv";
  if ($_POST["data_choice"]>1) {$dummy = explode('.',$dataset); $dataset=$dummy[0].'_ID'.$dummy[1];}
  break;
case 2:
  $go_file = "msigdb.txt";
  break;
default:
  $go_file = "human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv";
}

if ($help=='') {
$output = exec("python validate_gene_list.py $dataset $go_file");
$split_output = explode(",", $output);
//echo $dataset.','.$go_file; //for testing purposes
echo $split_output[0].','.$split_output[1].','.$split_output[2];
} else {echo "<ul id='error' style='color:red;'>$help</ul>";}
?>