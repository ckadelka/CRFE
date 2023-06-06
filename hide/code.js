$('.explain').hide();
        
$('input').mouseover(function(){ $(this).trigger('focus'); $(this).select();});

var choiceBox = $("input[name='alg_choice']:checked").val();    
changeChoiceBox(choiceBox);

$("input[name='alg_choice']").change(function() {
    changeChoiceBox($("input[name='alg_choice']:checked").val());
});

var paramLearningBox = $("input[name='param_learning_choice']");    
changeParamLearningBox(choiceBox);

$("input[name='param_learning_choice']").change(function() {
    changeParamLearningBox($("input[name='alg_choice']:checked").val());
});