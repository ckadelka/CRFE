$(document).ready(function() {
    $(".button").click(function() {
        document.getElementById("result").innerHTML = "Calculating..."; 
        window.scrollTo(0,99999);
        var formData = new FormData($("#myForm")[0]);
        
        $.ajax({
            type: "POST",
            url: "php_to_python_eval.php",
            data: formData,
            processData: false,
            contentType: false,
            success: function(data) {
                document.getElementById("result").innerHTML = data;
                window.scrollTo(0,99999);
                do_table_stuff();
            }
        });
        
        return false;
    });
    
    var last_chosen = 3;
    
    $('#file_gene_list').change(function(){
        var file = this.files[0];
        var name = file.name;
        var size = file.size;
        var type = file.type;
        if(name != '') { 
            var valid_extensions = /(.txt)$/i;   
            if(valid_extensions.test(name)) {
                var formData = new FormData($("#myForm")[0]);
                $.ajax({
                    type: "POST",
                    url: "php_to_python_validate_gene_list.php",
                    data: formData,
                    processData: false,
                    contentType: false,
                    success: function(data) {
                        var split_data=data.split(",");
                        if (parseInt(split_data[2])>0) {
                            alert("File accepted.\nOf the "+split_data[0]+" genes in the gene file, "+split_data[2]+" were found in the ontology." );
                            calculate_perturbed_genes(last_chosen);
                        } else {
                            alert("File accepted.\nHowever, none of the genes in the gene file are also in the ontology. Check that the same name space is used in ontology and gene list.");
                        }
                    }
                });
                
            }
            else {
                alert('Invalid File. Only txt files allowed.');
                document.getElementById('file_gene_list').value='';
            }
        } 
        
        return false;
    });


    $('.explain').hide();
    $("#textarea").hide();
    
    change_onto_explanation();
    change_data_kind_explanation();
        
    $('input').mouseover(function(){ $(this).trigger('focus'); $(this).select();});
    
    var choiceBox = $("input[name='alg_choice']:checked").val();    
    changeChoiceBox(choiceBox);
    
    $("input[name='alg_choice']").change(function() {
        changeChoiceBox($("input[name='alg_choice']:checked").val());
    });
    
    var weightchoiceBox = $("input[name='weight_choice']:checked").val();    
    changeWeightChoiceBox(weightchoiceBox);
    
    $("input[name='weight_choice']").change(function() {
        changeWeightChoiceBox($("input[name='weight_choice']:checked").val());
    });        
        
    var paramLearningBox = $("input[name='param_learning_choice']");    
    changeParamLearningBox(choiceBox);
    
    $("input[name='param_learning_choice']").change(function() {
        changeParamLearningBox($("input[name='alg_choice']:checked").val());
    });
    
    $('#show_histogram').click(function(){
        show_histogram("hist_cc_web.py","Number of terms with a given number of gene annotations",'Genes Annotated to Term','Number of Occurences');

        //show_histogram("php_to_hist_cc.php","Number of terms with a given number of gene annotations",'Genes Annotated to Term','Number of Occurences');
        });
        
    $('#show_histogram2').click(function(){
        show_histogram("hist_mm_web.py","Number of terms with a given average expression level",'Average Expression Level','Number of Occurences');
        });
        
    $('#show_histogram3').click(function(){
        show_histogram("hist_gg_web.py","Number of genes with a given number of annotations",'Number of Annotations','Number of Occurences');
        });
        
    $('#show_histogram4').click(function(){
        //show_columnchart_old("php_to_hist_gg_google.php","Number of genes with a given number of annotations",'Number of Annotations','Number of Occurences');
        show_columnchart("hist_gg_web_google.py","Number of genes with a given number of annotations",'Number of Annotations','Number of Occurences');
        });
    
    $('.change_data').click(function(){
        if ($("input[name='data_choice']:checked").val()==0) {
            $("#textarea").show();
            //document.getElementById("result").innerHTML="AAA"+document.getElementById("file_gene_list").files[0].name;
            if (document.getElementById("file_gene_list").files[0].name != '') {calculate_perturbed_genes(last_chosen);}
        }
        else {
            $("#textarea").hide();
            calculate_perturbed_genes(last_chosen);
        }
    });
    
    $('.change_data_kind').click(function(){
        change_data_kind_explanation();
    });
    
    $('.change_onto').click(function(){
        calculate_perturbed_genes(last_chosen);
        change_onto_explanation();
    });
    
    $('.change_data_kind').click(function(){
        calculate_perturbed_genes(last_chosen);   
    });
    
    $("input[name='cutoff']").change(function() {
        calculate_perturbed_genes(last_chosen);   
    });
    
    $("input[name='top_cutoff']").change(function() {
        calculate_perturbed_genes(last_chosen);   
    });                    
                                                                
    $("input[name='min_level']").change(function() {
        calculate_perturbed_genes(1);
        last_chosen = 1;
    });
    
    
    $("input[name='nr_perturbed']").change(function() {
        calculate_perturbed_genes(2);
        last_chosen = 2;
    });
    
    $("input[name='proportion_perturbed']").change(function() {
        calculate_perturbed_genes(3);
        last_chosen = 3;
    });
});

function do_table_stuff(){
    zebraRows('table.sortable tbody tr:odd td', 'odd');
	
    $('table.sortable tbody tr').hover(function(){
        $(this).find('td').addClass('hovered');},
    function(){
        $(this).find('td').removeClass('hovered');
    });
   	
    //default each row to visible
    $('table.sortable tbody tr').addClass('visible');
    
    //overrides CSS display:none property
    //so only users w/ JS will see the
    //filter box
    $('#search').show();
    
    $('#filter').keyup(function(event) {
	//if esc is pressed or nothing is entered
        if (event.keyCode == 27 || $(this).val() == '') {
            //if esc is pressed we want to clear the value of search box
            $(this).val('');
            //each row should be visible because if nothing is entered then all rows are matched.
            $('table.sortable tbody tr').removeClass('visible').show().addClass('visible');
        }

        //if there is text, lets filter
	else {
            filter('table.sortable tbody tr', $(this).val());
        }

        //reapply zebra rows
        $('.visible td').removeClass('odd');
        zebraRows('.visible:even td', 'odd');
    });
	
    //grab all header rows
    $('table.sortable thead th').each(function(column) {
        $(this).addClass('sortable').click(function(){
            var findSortKey = function(cell) {
                return cell.find('.sort-key').text().toUpperCase() + ' ' + cell.text().toUpperCase();
            };
						
            var sortDirection = $(this).is('.sorted-asc') ? -1 : 1;					
            //step back up the tree and get the rows with data for sorting
            var rows = $(this).parent().parent().parent().find('tbody tr').get();
						
            //loop through all the rows and find 
            $.each(rows, function(index, row) {
                row.sortKey = findSortKey($(row).children('td').eq(column));
            });
				
            //compare and sort the rows alphabetically
            rows.sort(function(a, b) {
                //if(! isNaN (a.sortKey - 0) && ! isNaN (b.sortKey - 0) && isInteger(a.sortKey) && isInteger(b.sortKey)) return (parseInt(a.sortKey) - parseInt(b.sortKey))*sortDirection;
                if(! isNaN (a.sortKey - 0) && ! isNaN (b.sortKey - 0)) return (parseFloat(a.sortKey) - parseFloat(b.sortKey))*sortDirection;
                if (a.sortKey < b.sortKey) return -sortDirection;
                if (a.sortKey > b.sortKey) return sortDirection;
                return 0;
            });
				
            //add the rows in the correct order to the bottom of the table
            $.each(rows, function(index, row) {
                $('table.sortable tbody').append(row);
                row.sortKey = null;
            });
				
            //identify the column sort order
            $('th').removeClass('sorted-asc sorted-desc');
            var sortHead = $('th').filter(':nth-child(' + (column + 1) + ')');
            sortDirection == 1 ? sortHead.addClass('sorted-asc') : $sortHead.addClass('sorted-desc');
				
            //identify the column to be sorted by
            $('td').removeClass('sorted').filter(':nth-child(' + (column + 1) + ')').addClass('sorted');
				
            $('.visible td').removeClass('odd');
            zebraRows('.visible:even td', 'odd');
        });
    });
};


//used to apply alternating row styles
function zebraRows(selector, className)
{
    $(selector).removeClass(className)
	.addClass(className);
};

//filter results based on query
function filter(selector, query) {
	query	=	$.trim(query); //trim white space
  query = query.replace(/ /gi, '|'); //add OR for regex
  
  $(selector).each(function() {
    ($(this).text().search(new RegExp(query, "i")) < 0) ? $(this).hide().removeClass('visible') : $(this).show().addClass('visible');
  });
};

function changeChoiceBox(choiceBox) {
    $('.explain').hide();
    $('.mcmc_specific').hide();
    $('.mcmc_param_learning_specific').hide();
    $('#explain_param_learning').show();
    if (choiceBox == "Greedy") {
        $('#explain_Greedy').show();
    } else if (choiceBox == "MCMC") {
        $('#explain_MCMC').show();
        $('.mcmc_specific').show();
        if (document.getElementsByName("param_learning_choice")[0].checked == true) {$('.mcmc_param_learning_specific').show();}
    } else {
        jQuery.error = console.error;
        jquery.error("wrong choice");
        return 1;
    }
};

function changeWeightChoiceBox(weightchoiceBox) {
    if (weightchoiceBox == "Linear") {
        $('#explain_Linear').show();
        $('#explain_Expression').hide();
    } else if (weightchoiceBox == "Expression") {
        $('#explain_Linear').hide();
        $('#explain_Expression').show();
    } else {
        jQuery.error = console.error;
        jquery.error("wrong choice");
        return 1;
    }
};

function changeParamLearningBox(choiceBox) {
    $('.parameters').show();
    $('.mcmc_param_learning_specific').hide();
    
    if (document.getElementsByName("param_learning_choice")[0].checked == true) {
        $('.parameters').hide();
        if (choiceBox == "MCMC") {$('.mcmc_param_learning_specific').show();}
        
    }
};

function change_onto_explanation() {
    $(".explain_onto").hide();
    if ($("input[name='onto_choice']:checked").val()==1 || $("input[name='data_choice']:checked").val()==1) {
        $("#explain_GO").show();
    } 
    else if ($("input[name='onto_choice']:checked").val()==2) {$("#explain_msigdb").show();}
    else {$("#explain_own_onto").show();}
};
function change_data_kind_explanation() {
    $(".explain_data_kind").hide();
    if ($("input[name='data_kind_choice']:checked").val()==0) {$("#explain_normal").show();} 
    else if ($("input[name='data_kind_choice']:checked").val()==1) {$("#explain_inverted").show();}
    else {$("#explain_absolute").show();}
};

function isInteger(value) {
    if ((undefined === value) || (null === value)) {
        return false;
    }
    return value % 1 == 0;
};

function calculate_perturbed_genes(kind) {
    if (kind==1) {var min_lev = $("input[name='min_level']").val(); document.getElementsByName('nr_perturbed')[0].value="Calculating"; document.getElementsByName('proportion_perturbed')[0].value="Calculating";}
    else if (kind==2) {var min_lev = $("input[name='nr_perturbed']").val(); document.getElementsByName('min_level')[0].value="Calculating"; document.getElementsByName('proportion_perturbed')[0].value="Calculating";}
    else if (kind==3) {var min_lev = $("input[name='proportion_perturbed']").val(); document.getElementsByName('min_level')[0].value="Calculating"; document.getElementsByName('nr_perturbed')[0].value="Calculating";}
    
    //var data_choice = $("input[name='data_choice']").val();
    var data_choice = $('input[name="data_choice"]:checked').val();
    var onto_choice = $('input[name="onto_choice"]:checked').val();

    //document.getElementsByName('s')[0].setAttribute("value",data_choice);
    if (data_choice==0) {var dataset = "uploaded_data/"+$("input[name='file_gene_list']").val().replace(/^.*[\\\/]/, '');}
    else if (data_choice==1) {var dataset = "bkz_gene_list.csv";}
    else if (data_choice==2) {var dataset = "3-cell-hep-d12_versus_2-cell-hep-d12.txt";}
    else if (data_choice==3) {var dataset = "3-cell-hep-d12_versus_CS-d12.txt";}
    else if (data_choice==4) {var dataset = "3-cell-hep-d12_versus_HM-d12.txt";}
    else if (data_choice==5) {var dataset = "2-cell-hep-d12_versus_CS-d12.txt";}
    else if (data_choice==6) {var dataset = "2-cell-hep-d12_versus_HM-d12.txt";}
    else if (data_choice==7) {var dataset = "CS-d12_versus_HM-d12.txt";}
    
    if (onto_choice==1) {var var_go_file = "human-taxon-id-gene2go-with-annotation-status-closed-bioprocess-only.csv";}
    else if (onto_choice==2) {var var_go_file = "msigdb.txt";}
    
    $.post("php_to_python.php", { gene_file: dataset, min_level: min_lev, cutoff:$("input[name='cutoff']").val(), top_cutoff:$("input[name='top_cutoff']").val(), k: kind, go_file: var_go_file, kind_of_list:$('input[name="data_kind_choice"]:checked').val()}, function(data){
        var split_data=data.split(",");
        document.getElementsByName('min_level')[0].value=split_data[0];
        document.getElementsByName('nr_perturbed')[0].value=split_data[1];
        document.getElementsByName('proportion_perturbed')[0].value=split_data[2];
    });

};

function show_histogram(pythonurl,title,xlabel,ylabel) {
   var to_plot = [];
   var options = {
        bars: { show: true, align: "left" },
        xaxis: {show : true, axisLabel: xlabel},
        yaxis: {show : true, axisLabel: ylabel}
    };
    
    function onDataReceived(series) {
        to_plot.push(series);
        document.getElementById("title_inModal").innerHTML = title;
        options.bars["barWidth"] = series["data"][1][0] - series["data"][0][0];
        if (series["data"][1][0] - series["data"][0][0] == 1) {
            options.bars["align"] = "center";}
        $.plot($("#inModal"), to_plot,options);
        $("#myModal").reveal($(this).data());
    }
         
    $.post("php_to_plot.php", { data_choice: $('input[name="data_choice"]:checked').val(), cutoff:$("input[name='cutoff']").val(), top_cutoff:$("input[name='top_cutoff']").val(), onto_choice: $('input[name="onto_choice"]:checked').val(), kind_of_list:$('input[name="data_kind_choice"]:checked').val(), python: pythonurl}, function(data) {onDataReceived(data);}, 'json');
};

function show_columnchart(pythonurl,str_title,xlabel,ylabel) {    
    var jsonData = $.ajax({
        type: "POST",
        url: "php_to_plot.php",
        data: {data_choice: $('input[name="data_choice"]:checked').val(), cutoff:$("input[name='cutoff']").val(), top_cutoff:$("input[name='top_cutoff']").val(), onto_choice: $('input[name="onto_choice"]:checked').val(), kind_of_list:$('input[name="data_kind_choice"]:checked').val(), python: pythonurl},
        dataType:"json",
        async: false
        }).responseText;
         
    // Create our data table out of JSON data loaded from server.
    var data = new google.visualization.DataTable(jsonData);
    
    var options = {
        title: str_title,
        hAxis: {title: xlabel},
        vAxis: {title: ylabel},
        tooltip:{trigger: 'hover'}};
    
    var chart = new google.visualization.ColumnChart(document.getElementById('inModal'));
    chart.draw(data, options);
    document.getElementById("title_inModal").innerHTML = "";
    $("#myModal").reveal($(this).data());
};