// Set the parameters
let params = {};
params.count = {
  system: 'hilo',
  tc: 0, // True Count
  decks: 3 // Decks remaining in shoe
}

// Construct Jackfish
let jackfish = new Jackfish(params, callback);

// Button onclick function
function setTrueCount(tc) {
  params.count.tc = tc;
  jackfish.setParams(params, true, callback);
}

// When jackfish is done computing, call this function
function callback(data) {
  let tableText = '';
  data.table.forEach(row => {
    row.forEach(cell => {
      tableText += cell.action;
    });
    tableText += '<br />';
  });
  tableText += '<br />';
  setHTML(tableText, `Edge: ${data.edge}<br />`, `Insurance: ${data.insurance}`);
}

// Set HTML in #output
function setHTML(...texts) {
  document.getElementById('output').innerHTML = '';
  texts.forEach((text) => {
    document.getElementById('output').innerHTML += text;
  })
}
