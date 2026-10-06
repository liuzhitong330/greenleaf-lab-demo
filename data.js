(function () {
  "use strict";
  window.GREENLEAF_DATA_READY = fetch("data/data.json?v=20261006")
    .then(function (response) {
      if (!response.ok) throw new Error("Count data could not be loaded (" + response.status + ").");
      return response.json();
    });
}());
