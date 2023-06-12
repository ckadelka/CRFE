<?php


$help = "";
$kind_of_list = $_POST["data_kind_choice"];
if ($_POST["data_choice"]==1 and $_POST["onto_choice"]==2) $help .= "<li>MSigDB cannot be used for HDF data yet. Please choose GO - Biological Processes.</li>";
$nr_categories = $_POST["nr_categories"];
if ((preg_match('/^\d+$/',$nr_categories)) == false || $nr_categories<0) $help .= "<li>Number of categories must be a nonnegative integer.</li>";
if ($_POST["alg_choice"] == "MCMC") $MCMC=1;
else $MCMC=0;
if ($_POST["weight_choice"] == "Linear") $weight_kind=1;
else $weight_kind=0;
$burnout= $_POST["burnout"];
if ((preg_match('/^\d+$/',$burnout)) == false || $burnout<0) $help .= "<li>A nonnegative integer must be chosen for the burnout period.</li>";

$MCMC_steps= $_POST["MCMC_steps"];
if ((preg_match('/^\d+$/',$MCMC_steps)) == false || $MCMC_steps<=0) $help .= "<li>A positive integer must be chosen for the number of MCMC iterations.</li>";
$alpha_beta_top= $_POST["alpha_beta_top"];
if (is_numeric($alpha_beta_top) == false || $alpha_beta_top<=0 || $alpha_beta_top>1) $help .= "<li>A value in the interval (0,1] must be chosen as maximal value for alpha and beta.</li>";
$percentage_parameter_change= $_POST["percentage_parameter_change"];
if (is_numeric($percentage_parameter_change) == false || $percentage_parameter_change<=0 || $percentage_parameter_change>1) $help .= "<li>A value in the interval (0,1] must be chosen as the proportion of time for parameter change.</li>";
if ($_POST["param_learning_choice"]==false) $param_learning=0;
else $param_learning=1;
$alpha = $_POST['alpha'];
if (is_numeric($alpha) == false || $alpha<=0 || $alpha>=1) $help .= "<li>A value in the interval (0,1) must be chosen for the false positive rate alpha.</li>";
$beta = $_POST["beta"];
if (is_numeric($beta) == false || $beta<=0 || $beta>=1) $help .= "<li>A value in the interval (0,1) must be chosen for the false negative rate beta.</li>";
$prob = $_POST["prob"];
if (is_numeric($prob) == false || $prob<=0 || $prob>=0.5) $help .= "<li>A value in the interval (0,0.5) must be chosen for the penalization parameter prob.</li>";
$cutoff = $_POST["cutoff"];
if ((preg_match('/^\d+$/',$cutoff)) == false || $cutoff<=0) $help .= "<li>A positive integer must be chosen for the minimal number of genes annotated by each term.</li>";
$top_cutoff = $_POST["top_cutoff"];
if ((preg_match('/^\d+$/',$top_cutoff)) == false || $top_cutoff<0) $help .= "<li>A positive integer must be chosen for the maximal number of genes annotated by each term.</li>";
else if ((preg_match('/^\d+$/',$cutoff)) == true && $top_cutoff > 0 && $top_cutoff < $cutoff) $help .= "<li>The maximal number of genes annotated by each term cannot be smaller than the minimal number of genes annotated by each term.</li>";
$belief = $_POST["belief"];
if (is_numeric($belief) == false || $belief<=0) $help .= "<li>A positive value must be chosen for the belief parameter. As long as you trust the experimental data at all, the belief parameter should be greater than 1.</li>";
$min_level = $_POST["min_level"];
if (is_numeric($min_level) == false) $help .= "<li>A numeric value must be chosen for the threshold between perturbed and unperturbed genes.</li>";
$nr_perturbed = $_POST["nr_perturbed"];
if ((preg_match('/^\d+$/',$nr_perturbed)) == false || $nr_perturbed<0) $help .= "<li>The number of perturbed genes must be a nonnegative integer.</li>";
$proportion_perturbed = $_POST["proportion_perturbed"];
if (is_numeric($proportion_perturbed) == false || $proportion_perturbed<0 || $proportion_perturbed>1) $help .= "<li>The ratio of perturbed to total genes must be any value in the interval [0,1].</li>";

$extension=".txt";

switch ($_POST["data_choice"])
{
case 0:
    $allowedExts = array("txt");
    $temp = explode(".", $_FILES["file_gene_list"]["name"]);
    $file_extension = end($temp);
    $file=$_FILES["file_gene_list"];
    $len_temp = sizeof($temp[0]);
    //
    //echo "<p>FILE: $file $len_temp</p>";
    //echo "<p>Extension: $extension Allowed: $allowedExts TEMP: $temp $temp[0] $temp[1]</p>";
    
    if (($_FILES["file_gene_list"]["type"] == "text/plain")
    && ($_FILES["file_gene_list"]["size"] < 2000000)
    && in_array($file_extension, $allowedExts)) {
    if ($_FILES["file_gene_list"]["error"] > 0) {
        echo "Error: " . $_FILES["file_gene_list"]["error"] . "<br>";
        $help .="<li>Error with the data file: ".$_FILES["file_gene_list"]["error"]."</li>";
    } else {
        echo "Upload: " . $_FILES["file_gene_list"]["name"] . "<br>";
        echo "Type: " . $_FILES["file_gene_list"]["type"] . "<br>";
        echo "Size: " . ($_FILES["file_gene_list"]["size"] / 1024) . " kB<br>";
        echo "Stored in: " . $_FILES["file_gene_list"]["tmp_name"];
        //$uploaddir = '/var/www/uploads/';
        $uploadfile = basename($_FILES['file_gene_list']['name']);
        //echo "$uploadfile";
        
        echo '<pre>';
        if (move_uploaded_file($_FILES['file_gene_list']['tmp_name'], "uploaded_data/".$uploadfile)) {
            echo "File is valid, and was successfully uploaded.\n";
            $dataset = "uploaded_data/".$temp[0];
        } else {
            echo "Possible file upload attack!\n";
        }
        echo '</pre>';
    }
    } else {
    echo "<p>Invalid file</p>";
    $help .= "<li>Invalid file for data. Check that you chose a txt file.</li>";
    }
    break;
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
  if ($_POST["data_choice"]>1) {$dummy = explode('.',$dataset); $dataset=$dummy[0].'_ID'.$dummy[1];}
  break;
case 2:
  $go_file = "msigdb.txt";
  break;
default:
  $go_file = "human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv";
}
$dataset = $dataset.$extension; 

$outstr = ($_POST["data_choice"] > 0) ? (($_POST["data_choice"] > 1) ? "LIVER" : "HDF") : "OWN";

/*
if ($help!="") {echo "<ul id='error' style='color:red;'>$help</ul>";}
else {echo $burnout.','.$MCMC.','.$weight_kind.','.$outstr.','.$go_file.','.$dataset.','.$weight_kind.','.$param_learning.','.$alpha;
}*/

//echo "AAAAAAAA -L$nr_categories -M$MCMC -P$param_learning -a$alpha -b$beta -p$prob -l$belief";

if ($help!="") {echo "<ul id='error' style='color:red;'>$help</ul>";}
else {
    $command = "python GeneSetRankedEnrichment32e_web.py";
    $command .= " -L$nr_categories -M$MCMC -P$param_learning -a$alpha -b$beta -p$prob -l$belief -m$proportion_perturbed -c$cutoff -z$top_cutoff -e$dataset -d$outstr -B$burnout -S$MCMC_steps -X$alpha_beta_top -E$kind_of_list";
    
    echo "$command";
    
    $pid = popen( $command,"r");
        
        while( !feof( $pid ) )
        {
            echo fread($pid, 1000);
            flush();
            ob_flush();
            //echo "<script>window.scrollTo(0,99999);</script>";
            usleep(100000);
        }
    pclose($pid);
}
?>