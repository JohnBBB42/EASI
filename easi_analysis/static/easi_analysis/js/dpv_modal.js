// Get elements
var analysisSelect = document.getElementById("id_selected_analysis");
var modal = document.getElementById("plotModal");
var closeBtn = document.querySelector("#plotModal .close");
var selectAllBox = document.getElementById("selectAll");
var applyBtn = document.getElementById("applySelection");

// Show modal when the specific analysis option is selected
analysisSelect.addEventListener("change", function() {
  if (
    this.value === "Plot Data Array with Corrected Baseline" ||
    this.value === "Analyze Peak Currents" ||
    this.value === "Observed vs Expected Concentration"
  )

    modal.style.display = "flex";  // show modal (using flex for centering)
  } else {
    // If switching to a different analysis, ensure modal is hidden
    modal.style.display = "none";
    // (Optional) clear any checked boxes if analysis changed:
    selectAllBox.checked = false;
    document.querySelectorAll('input[name="selected_plots"]').forEach(cb => {
      cb.checked = false;
    });
  }
});

// Close modal when clicking the close X or the OK button
closeBtn.addEventListener("click", function() {
  modal.style.display = "none";
});
applyBtn.addEventListener("click", function() {
  modal.style.display = "none";
});

// Close modal if clicking outside the modal content
window.addEventListener("click", function(event) {
  if (event.target === modal) {
    modal.style.display = "none";
  }
});

// "Select All" functionality: toggle all plot checkboxes
selectAllBox.addEventListener("change", function() {
  var plotCheckboxes = document.querySelectorAll('input[name="selected_plots"]');
  plotCheckboxes.forEach(cb => {
    cb.checked = selectAllBox.checked;
  });
});
