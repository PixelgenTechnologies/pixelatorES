// prevent nav-pills jump
document.addEventListener("DOMContentLoaded", function() {
  document.querySelectorAll('.nav-pills a').forEach(function(a) {
    a.addEventListener('click', function(e) {
      e.preventDefault(); // prevent the page from jumping
    });
  });
});