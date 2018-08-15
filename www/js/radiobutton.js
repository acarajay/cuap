$(document).ready(function() {

	// Funtion to get checked radio's value.
	$('#radio_value').click(function() {
		$('#result').empty();
		var value = $("form input[type='radio']:checked").val();
		if($("form input[type='radio']").is(':checked')) {
			//$('#result').append(value);
			Shiny.onInputChange("radioval",value);
		}else{
			//value = "zero";
			//Shiny.onInputChange("radioval",value);
			alert(" Please click any option! ");
			location.reload();
		}
	});

	// Get value Onchange radio function.
	$('input:radio').change(function(){
		var value = $("form input[type='radio']:checked").val();
		Shiny.onInputChange("radioval",value);
		//alert("Value of Changed Radio is : " +value);
	});

	// Funtion to reset or clear selection.
	$('#form_reset').click(function() {
		//alert("Form reset");
		$('#result').empty();
		$("input:radio").attr("checked", false);
		location.reload();
	});

	$('#inputfile').filestyle({ 
		buttonName : 'btn-success',		 
		buttonText : ' Open'	 
	});
});