# Jackfish
A JavaScript Blackjack engine that determines:
 * Player/House edge
 * Perfect strategy

for standard blackjack given:
 * Rule variants
 * Number of decks left in shoe
 * Counting Systems
 * True Counts

Jackfish also includes a simulator that generates player edge and outcome frequencies, given:
 * A betting system based on the true count
 * Deck penetration

## License
This software is available under the MIT license.

## Installation
Clone this repository and add `Jackfish.js` and `JackfishWorker.js` to your project directory.
Add the following line of code to your html file:

```xml
<script type='text/javascript' src='Jackfish.js'></script>
```

The worker file will be included automatically by `Jackfish.js`

## Example
```JavaScript
let jackfish = new Jackfish({
  surrender: 'late',
  count: {
    system: 'hilo',
    count: 12,
    tc: 4, // True Count
    decks: 3
  }
});

jackfish.doAll(() => {
  jackfish.getEdge();                 // 0.016476440033674107
  jackfish.takeInsurance();           // true
  jackfish.getHit(16, 9);             // -0.5538886717476285
  jackfish.getStand(16, 9);           // -0.5537391635562128
  jackfish.getTable(16, 9).action;    // "S"
  jackfish.getTable(16, 9).surrender; // true
});
```

## Contributions
Contributions are open and welcome! There are many more rule variants and game variants (eg. Free Bet Blackjack) that Jackfish has not yet implemented. If your local casino uses a variant that Jackfish does not currently support, feel free to create a pull request!
