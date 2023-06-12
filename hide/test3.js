//$("#inModal").load("",function(){
    var data_choice = $('input[name="data_choice"]:checked').val();
    var onto_choice = $('input[name="onto_choice"]:checked').val();

    //document.getElementsByName('s')[0].setAttribute("value",data_choice);
    if (data_choice==1) {var dataset = "bkz_gene_list.csv";}
    else if (data_choice==2) {var dataset = "3-cell-hep-d12_versus_2-cell-hep-d12.txt";}
    else if (data_choice==3) {var dataset = "3-cell-hep-d12_versus_CS-d12.txt";}
    else if (data_choice==4) {var dataset = "3-cell-hep-d12_versus_HM-d12.txt";}
    else if (data_choice==5) {var dataset = "2-cell-hep-d12_versus_CS-d12.txt";}
    else if (data_choice==6) {var dataset = "2-cell-hep-d12_versus_HM-d12.txt";}
    else if (data_choice==7) {var dataset = "CS-d12_versus_HM-d12.txt";}
    
    if (onto_choice==1) {var var_go_file = "human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv";}
    else if (onto_choice==2) {var var_go_file = "msigdb.txt";}
    
    var jsonData = $.ajax({
        url: "php_to_hist_gg_google.php",
        data: {gene_file: dataset, cutoff:$("input[name='cutoff']").val(), top_cutoff:$("input[name='top_cutoff']").val(), go_file: var_go_file, kind_of_list:$('input[name="data_kind_choice"]:checked').val()},
        dataType:"json",
        async: false
        }).responseText;
  
    document.getElementById('result').innerHTML=jsonData;
        
    // Create our data table out of JSON data loaded from server.
    var data = new google.visualization.DataTable(jsonData);
    
    var options = {
        title: 'Company Performance',
        hAxis: {title: 'Questions', titleTextStyle: {color: 'red'}},
        vAxis: {title: '1 = POOR, 5 = EXCELLENT', titleTextStyle: {color: '#FF0000'}, maxValue:'5', minValue:'1'},
        tooltip:{trigger: 'hover'}};
    
    var chart = new google.visualization.ColumnChart(document.getElementById('inModal'));
    chart.draw(data, options);
    $("#myModal").reveal($(this).data());
//});