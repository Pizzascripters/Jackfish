function Jackfish(params, callback) {
  /*-- Constants --*/

  const CARD_ODDS = 1/13; // The chance that we draw a non-ten card
  const TEN_ODDS = 4/13;
  const ACE_HIGH_RATIO = 1/5; // Aces : high cards

  // Special hand cases
  const BLACKJACK = -1;
  const BUST = -2;
  const DEALER_TEN = -3;
  const DEALER_ACE = -4;

  const VALUE_MASK = 0x1f;
  const SOFT = 0x20;
  const PAIR = 0x40;
  const ACE = 11 | SOFT;
  const ACE_PAIR = ACE | PAIR;

  // A list of all cards
  const CARD_STATES = (() => {
    let i = 2;
    let a = fillArray(() => i++, 9);
    a.push(ACE);
    return a;
  })();
  // A list of possible hands the player or the dealer can have
  const HAND_STATES = (() => {
    let a = [];
    loop(2, 22, v => a.push(createHand(v, false))); // Hard hands
    loop(11, 21, v => a.push(createHand(v, true))); // Soft hands
    a.push(0, BLACKJACK, BUST, DEALER_TEN, DEALER_ACE);
    return a;
  })();
  // The order we should determine the best moves
  const HAND_ORDER = (() => {
    let a = [];
    a.push(BUST, BLACKJACK);
    rLoop(21, 11, v => a.push(createHand(v))); // Hard hands 12+
    rLoop(20, 10, v => a.push(createHand(v, true))); // Soft hands
    rLoop(11, 1, v => a.push(createHand(v))); // Hard hands <12
    a.push(0);
    loop(2, 11, v => a.push(createHand(v, false, true))); // Pairs
    a.push(createHand(ACE, false, true));
    return a;
  })();
  // The player hands recorded by the table
  const TABLE_HANDS = (() => {
    let a = [];
    loop(5, 22, v => a.push(createHand(v, false))); // Hard hands
    loop(12, 21, v => a.push(createHand(v, true))); // Soft hands
    loop(2, 11, v => a.push(createHand(v, false, true))); // Pairs
    a.push(createHand(ACE, false, true));
    return a;
  })();

  // All possible hands the dealer can start with
  const DEALER_STATES = HAND_STATES.filter(hand =>
    (hand > 0 && hand < 10) ||
    hand === DEALER_TEN ||
    hand === DEALER_ACE
  );

  const INDEX_BLACKJACK = HAND_STATES.indexOf(BLACKJACK);
  const INDEX_21 = HAND_STATES.indexOf(21);

  const P = progressionMatrix();

  /*-- Counting Constants --*/

  const HILO = [1, 1, 1, 1, 1, 0, 0, 0, -1, -1];
  const KO = [1, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1];
  const OMEGA2 = [1, 1, 2, 2, 2, 1, 0, -1, -2, 0];
  const WONGHALVES = [.5, 1, 1, 1.5, 1, .5, 0, -.5, -1, -1];
  const USTONAPC = [1, 2, 2, 3, 2, 2, 1, -1, -3, 0];
  const SYSTEMS = {
    'hilo': HILO,
    'ko': KO,
    'omega2': OMEGA2,
    'wonghalves': WONGHALVES,
    'ustonapc': USTONAPC
  }

  /*-- Private variables --*/

  let comp, // 10d array representing the composition of the shoe
      table, // Matrix indicating the best move given dealer up card and player hand
      matricesMade = false;

  // Return matrices. Specifies return given player's hand and dealer's card under perfect play
  let rM = zeroes([HAND_STATES.length, DEALER_STATES.length]); // Return without doubling
  let rdM = zeroes([HAND_STATES.length, DEALER_STATES.length]); // Return with doubling
  let rsM = zeroes([HAND_STATES.length, DEALER_STATES.length]); // Return with surrendering
  let standM = zeroes([HAND_STATES.length, DEALER_STATES.length]);
  let hitM = zeroes([HAND_STATES.length, DEALER_STATES.length]);
  let doubleM = zeroes([HAND_STATES.length, DEALER_STATES.length]);
  let splitM = zeroes([DEALER_STATES.length, DEALER_STATES.length]);

  /*-- Public functions --*/

  this.doAll = (callback) => {
    let all = {};
    let valuesLeft = 4;

    this.makeMatrices(subCallback.bind(null, 'matrices'));
    this.makeTable(subCallback.bind(null, 'table'));
    this.takeInsurance(subCallback.bind(null, 'insurance'));
    this.getEdge(subCallback.bind(null, 'edge'));

    function subCallback(key, value) {
      all[key] = value;
      if(--valuesLeft <= 0) {
        all.bjOdds = DEALER_STATES.map(getBJOdds.bind(this, comp));
        all.matrices.rM = rM;
        all.matrices.rdM = rdM;
        all.matrices.rsM = rsM;
        all.matrices.hitM = hitM;
        all.matrices.splitM = splitM;
        if(callback) {
          callback(all);
        }
      }
    }
  }

  this.setParams = (callback, params_) => {
    params = params_;
    if(params.count.system !== 'none') {
      params.count.indices = SYSTEMS[params.count.system];
    }
    this.comp = comp = deckComp(params.count);
    if(callback) {
      callback(params);
    }
  }

  // Make endM, standM, and doubleM
  this.makeMatrices = (callback, comp) => {
    if(comp === undefined) {
      comp = this.comp;
    }

    // Calculate dealer's odds to reach each endstate
    endM = endMatrix(comp, params.soft17, params.count.cards);

    // Calculate player's return by standing
    standM = standReturns(endM);

    // Set returns on 21, bust, and blackjack
    DEALER_STATES.forEach((dealer, j) => {
      // Determine return on blackjack
      let bjReturn = params.blackjack;
      if(!params.peek && (dealer === DEALER_TEN || dealer === DEALER_ACE)) {
        // No peek; we can push if dealer blackjack
        bjReturn *= 1 - bjOdds;
      }

      // Set the returns. [hand, return]
      [[21], [BUST, -1], [BLACKJACK, bjReturn]].forEach((a) => {
        let i = HAND_STATES.indexOf(a[0]);
        if(a[1] === undefined) {
          a[1] = standM[i][j];
        }
        rsM[i][j] = rdM[i][j] = rM[i][j] = a[1];
      });
    });

    // Calculate player's return by doubling
    doubleM = doubleReturns(comp, standM);

    matricesMade = true;
    if(callback) {
      callback({endM, standM, doubleM});
    }
    return {endM, standM, doubleM};
  }

  this.makeTable = (callback, comp) => {
    if(comp === undefined) {
      comp = this.comp;
    }
    table = [];

    // Search the cache to see if this table was already made
    let cache = searchTableCache(params);
    if(loadCache(searchTableCache(params))) {
      // Table was previously cached
      if(callback) {
        callback(table);
      }
      return table;
    }

    this.makeMatrices(null, comp);

    HAND_ORDER.forEach(player => {
      if(player < 0) return;
      let pair = player & PAIR;
      let i = HAND_STATES.indexOf(player);
      let k = TABLE_HANDS.indexOf(player);
      if(k !== -1) {
        table[k] = [];
      }
      DEALER_STATES.forEach((dealer, j) => {
        function doCell() {
          let move = bestMove(comp, player & 0x3f, dealer, pair);
          if(!pair) {
            rsM[i][j] = move.expected;
            rdM[i][j] = move.expectedNS;
            rM[i][j] = Math.max(standM[i][j], hitM[i][j]);
          }
          if(k !== -1) {
            table[k][j] = move;
          }
        }
        doCell();
      });
    });

    if(callback) {
      callback(table);
    }

    cacheTable(params, table);

    return table;
  }

  /*
   * Determines whether player should take insurance
   * cards: {count: Count, decks: Number, trueCount: Number}
   * Deck composition is derived from cards
   * If cards === undefined, comp is used
   * */
  this.takeInsurance = (callback, cards) => {
    // If cards !== undefined,
    let tempComp = comp;
    let decks = params.count.decks;
    if(cards !== undefined && cards.count !== undefined) {
      decks = cards.cards / 52;
      tempComp = countingComp({
        indices: params.count.indices,
        tc: cards.trueCount,
        decks
      });
    }
    if(params.count.system === 'hilo' || params.count.system === 'wonghalves') {
      /*
       * In systems where Aces and Tens are counted the same,
       * There's one more ten and one less ace than we think there is
       * Because the dealer has an Ace showing face up
       * */
      insurance = tempComp[CARD_STATES.indexOf(10)] + 1/(52 * decks) > 1/3;
    } else {
      insurance = tempComp[CARD_STATES.indexOf(10)] > 1/3;
    }
    if(callback) {
      callback(insurance);
    }
    return insurance;
  }

  // If j !== undefined, give the edge for dealer card j
  // Otherwise, give total edge
  this.getEdge = (cb, j) => {
    if(j === undefined) {
      let edge = 0;
      DEALER_STATES.forEach((_, j) => {
        edge += comp[j] * this.getEdge(null, j);
      });
      cb(edge);
      return edge;
    }

    // Get the vector state after 2 cards
    let shiftedComp = pullCard(comp, j, 52 * params.count.decks);
    let dealer = DEALER_STATES[j],
        state = squishMatrix(progressMState(
          initStateMatrix(0, shiftedComp),
          52 * params.count.decks, [], 2
        ));

    // Add blackjack to state
    state[INDEX_BLACKJACK] = state[INDEX_21]; // 21 after 2 cards is blackjack
    state[INDEX_21] = 0;

    // Remove pairs from state
    let rPair = []; // 10d vector of the return on each pair (no surrender)
    CARD_STATES.forEach((card, c) => {
      let i = HAND_STATES.indexOf(pairState(card));
      let pairOdds = shiftedComp[c] * shiftedComp[c];
      let rNS = table[TABLE_HANDS.indexOf(card | PAIR)][j].expectedNS;
      state[i] -= pairOdds; // 8,8 is a different hand than 6,10 so subtract from state
      rPair.push(pairOdds * rNS); // Store the return in rPair
    });

    // Total return is the state dotted with surrender returns plus pair returns
    let rTotal = dot(state, transpose(rsM)[j]) + vtotal(rPair);

    // Calculate return of insurance
    let insurance = 0;
    if(dealer === DEALER_ACE && this.takeInsurance()) {
      insurance = 1.5 * shiftedComp[CARD_STATES.indexOf(10)] - .5;
    }

    // Factor in dealer blackjack and return
    if(params.peek) {
      let bjOdds = getBJOdds(shiftedComp, CARD_STATES[j]);
      return rTotal * (1 - bjOdds) - bjOdds + insurance;
    } else {
      // If no peek, dealer blackjack was already factored into rTotal
      return rTotal + insurance;
    }
  }

  this.createSimulation = createSimulation;

  /*-- Monte Carlo Simulation --*/

  function createSimulation(options) {
    let cash, hands, shoes, // Assigned in reset
        totalHands, totalBet, totalDiff, frequencies, endRecord, mean, // Assigned in clear
        p, activeComp, tables; // Assigned in config
    let running = false,
        sessions = 0;

    let sim = {
      config: (options_) => {
        options = options_;
        if(options.player) {
          p = elemsToIndices(player);
        } else {
          p = false;
        }
        // Create card object
        cards = {
          cards: Math.round(52 * params.count.decks),
        };
        if(params.count.system !== 'none') {
          cards.count = 0;
          cards.trueCount = 0;
        }
        activeComp = copy(comp); // Real deck composition in game
        sim.clear();

        // Generate strategy tables to reference
        let tables = [];
        if(params.count.system !== 'none') {
          let originalTc = params.count.tc;
          for(let tc = -12; tc <= 12; tc++) {
            params.count.tc = tc;
            this.setParams(null, params);
            tables.push({
              trueCount: tc,
              table: this.makeTable()
            });
          }
          params.count.tc = originalTc;
          this.setParams(null, params);
        } else {
          tables = [{
            trueCount: 0,
            table: this.makeTable(null)
          }]
        }
      },

      run: (cb) => {
        running = true;
        let dealer = options.dealer === undefined ? undefined : DEALER_STATES.indexOf(options.dealer);

        while(
          cash > 0 &&
          hands < options.maxHands &&
          shoes < options.maxShoes &&
          (cash < options.maxCash || options.maxCash === Infinity)
        ) {
          let bet = options.bet;
          let maxTc = -Infinity;
          if(cards.trueCount !== undefined) {
            options.rules.forEach(rule => {
              if(cards.trueCount > maxTc && cards.trueCount > rule.value) {
                maxTc = rule.value;
                bet = rule.bet;
              }
            });
          }
          let r = playHand.bind(this)(activeComp, cards, {
            dealer,
            player: p,
            forceMove: options.forceMove,
            shuffledComp: comp,
            tables
          });
          hands++;
          totalHands++;
          totalBet += bet;
          totalDiff += r * bet;
          cash += r * bet;
          if(totalBet === 0) {
            mean = 0;
          } else {
            mean = totalDiff / totalBet;
          }
          if(bet > 0) {
            if(!frequencies[r]) {
              frequencies[r] = 1;
            } else {
              frequencies[r]++;
            }
          }
          if(cards.cards / 52 < options.yellow) {
            // Recopying comp simulates shuffle
            activeComp = shuffle(comp, cards);
            shoes++;
          }
        }

        let ending;
        if(cash <= 0) {
          ending = '$0';
        } else if(cash >= options.maxCash && cash !== Infinity) {
          ending = `$${options.maxCash}`;
        } else if(hands >= options.maxHands) {
          ending = `${options.maxHands} hands`;
        } else if(shoes >= options.maxShoes) {
          ending = `${options.maxShoes} shoes`;
        }
        if(ending && endRecord[ending]) {
          endRecord[ending]++;
        } else if(ending) {
          endRecord[ending] = 1;
        }

        setTimeout(() => {
          if(running) {
            sessions++;
            cb(sim.get());
            sim.reset();
            sim.run(cb);
          }
        }, 0);
      },

      stop: () => {
        running = false;
      },

      get: () => {
        return {
          cash,
          frequencies,
          endRecord,
          mean,
          totalHands
        }
      },

      getOptions: () => options,

      reset: () => {
        cash = options.cash;
        hands = 0;
        shoes = 0;
      },

      clear: () => {
        cash = options.cash;
        hands = 0;
        shoes = 0;
        frequencies = {};
        endRecord = {};
        mean = 0;
        totalHands = 0;
        totalBet = 0;
        totalDiff = 0;
      }
    }

    sim.config(options);

    // Converts an array of elements of CARD_STATES to indices
    function elemsToIndices(p) {
      if(!p) return p;
      let a = [];
      for(let i = 0; i < p.length; i++) {
        a[i] = CARD_STATES.indexOf(p[i]);
      }
      return a;
    }

    return sim;
  }

  // comp is a column vector representing the deck composition
  // cards are number of cards left in deck
  // options.dealer is the index of dealer card
  // options.player is a length 2 vector with the indices of the player cards
  // options.fixDealer prevents the dealer from hitting
  // options.forceMove dictates which move the player must play
  // options.noDouble
  let ten = 0;
  let ace = 0;
  let bj = 0;
  let total = 0;
  function playHand(comp, cards, options) {
    if(!options) {
      options = {};
    }

    let dealer,   // Dealer hand
        player,   // Player hand
        pair,     // Boolean: if player hand is a pair
        half;     // Half of the player value
    let double = false; // Whether or not we doubled

    // If using reference tables, find the appropriate table based on the count
    function getActiveTable() {
      if(options.tables && cards.trueCount) {
        let minDistance = Infinity;
        let closestTable = null;
        options.tables.forEach((table) => {
          let distance = Math.abs(table.trueCount - cards.trueCount);
          if(distance < minDistance) {
            minDistance = distance;
            closestTable = table;
          }
        });
        return closestTable.table;
      } else if(options.tables) {
        return options.tables[0].table;
      } else {
        return null;
      }
    }

    // Draw player cards
    let oldCount = cards.trueCount;
    let p = []; // List of player cards
    if(!options.player || options.player.length === 0) {
      p[0] = drawCard(comp, cards, options.shuffledComp);
      p[1] = drawCard(comp, cards, options.shuffledComp);
    } else {
      p = copy(options.player);
    }

    // Determine player hand
    player = HAND_STATES.indexOf(0); // Stand with hand value of zero
    for(let card of p) {
      player = HAND_STATES.indexOf(P[player][card]);
    }
    player = HAND_STATES[player];
    pair = (p.length === 2) && (p[0] === p[1]);

    // Draw random dealer card if not specified
    if(options.dealer === undefined || options.dealer === -1) {
      options.dealer = drawCard(comp, cards, options.shuffledComp);
    }
    dealer = DEALER_STATES[options.dealer];

    // Handle blackjack and pairs
    if(player === 21 && p.length === 2) {
      // Player blackjack
      player = BLACKJACK;
    } else if(pair && player === 44) {
      // Pair of aces
      half = 11;
      player = createHand(half, true, true);
    } else if(pair) {
      // Pair of anything but aces
      half = player / 2;
      player = createHand(half, false, true);
    }

    // Player takes insurance?
    let insurance = dealer === DEALER_ACE && this.takeInsurance(a=>0, cards) ? true : 0;

    if(player !== BLACKJACK && player !== BUST && player !== 21) {
      // Determine action
      let action;
      if(options.forceMove) {
        action = options.forceMove;
      } else {
        action = getTable(player, dealer, getActiveTable()).action;
      }

      // Handle no doubling
      // TODO: Remove this when safe
      if(options.noDouble && action === 'D') {
        action = 'H';
      } else if(options.noDouble && action === 'd') {
        action = 'S';
      }

      // Do said action
      if(action === 'P') {
        if(!options.afterSplit) {
          dealer = simDealer(P, comp, cards, dealer, options.shuffledComp);
          if(insurance && dealer === BLACKJACK) {
            insurance = 1;
          } else if(insurance) {
            insurance = -.5;
          }
        } else if(insurance) {
          // Can't take insurance again after a split
          insurance = 0;
        }
        // Fix dealer card because all splits happen on the same table
        let options_ = {
          player: [p[0], drawCard(comp, cards, options.shuffledComp)],
          dealer: options.dealer,
          fixDealer: dealer,
          afterSplit: true,
          noDouble: !params.split.double,
          shuffledComp: options.shuffledComp
        };
        let r1 = playHand.bind(this)(comp, cards, options_);
        options.player = [p[1], drawCard(comp, cards, options.shuffledComp)]
        let r2 = playHand.bind(this)(comp, cards, options_);
        if(r1 === params.blackjack) r1 = 1; // Blackjack after split isn't natural
        if(r2 === params.blackjack) r2 = 1;
        return r1 + r2 + insurance;
      }
      if(pair) {
        player = P[HAND_STATES.indexOf(player & 0x3f)][HAND_STATES.indexOf(half)];
      }
      if(action === 'D' || action === 'd') {
        double = true;
        player = P[HAND_STATES.indexOf(player & 0x3f)][drawCard(comp, cards, options.shuffledComp)];
      } else {
        // Hit until table says not to
        while(action === 'H' || action === 'D') {
          player = P[HAND_STATES.indexOf(player & 0x3f)][drawCard(comp, cards, options.shuffledComp)];
          if(player !== BUST) {
            action = getTable(player, dealer, getActiveTable()).action;
          } else {
            break;
          }
          if(params.double.anytime && (action === 'D' || action === 'd')) {
            double = true;
            player = P[HAND_STATES.indexOf(player & 0x3f)][drawCard(comp, cards, options.shuffledComp)];
            break;
          }
        }
      }
    }

    dealer = simDealer(P, comp, cards, dealer, options.shuffledComp);

    // If insurance was taken, add or take away money
    if(insurance && dealer === BLACKJACK) {
      insurance = 1;
    } else if(insurance) {
      insurance = -.5;
    }

    if(player === BUST && double) {
      // Twice the losses for twice the bet
      return -2 + insurance;
    } else if(player === BUST) {
      return -1 + insurance;
    } else if(options.fixDealer) {
      dealer = options.fixDealer;
    }

    if(player === BLACKJACK && dealer !== BLACKJACK) {
      return params.blackjack + insurance; // Player blackjack
    } else if(params.peek && double && dealer === BLACKJACK) {
      return -1 + insurance;
    }

    let pv = getHandDetails(player)[0]; // Player value
    let dv = getHandDetails(dealer)[0]; // Dealer value
    if(pv > dv && double) {
      return 2 + insurance;
    } else if(pv > dv) {
      return 1 + insurance;
    } else if(pv === dv) {
      return insurance;
    } else if(!double) {
      return -1 + insurance;
    } else {
      return -2 + insurance;
    }
  }

  function drawCard(comp, cards, shuffledComp) {
    let tempComp = comp;
    if(shuffledComp && cards.cards <= 0) {
      tempComp = shuffle(shuffledComp, cards);
    }
    let c = pickFromArray(tempComp);
    let newComp = pullCard(tempComp, c, cards.cards);
    for(let i = 0; i < comp.length; i++) {
      comp[i] = newComp[i];
    }
    if(cards.count !== undefined) {
      cards.count += params.count.indices[c];
      cards.trueCount = 52 * cards.count / cards.cards;
    }
    cards.cards--;
    return c;
  }

  // Run dealer until it hits an endstate or bust
  function simDealer(P, comp, cards, dealer, shuffledComp) {
    while(
      (dealer & VALUE_MASK) < 17 ||
      dealer === DEALER_TEN ||
      dealer === DEALER_ACE ||
      (params.soft17 && dealer === 49)
    ) {
      let nextCard = drawCard(comp, cards, shuffledComp);
      let oldDealer = dealer;
      dealer = P[HAND_STATES.indexOf(dealer)][nextCard];
    }
    return dealer;
  }

  /*-- Private functions --*/

  let tableCache = [];
  function searchTableCache(params) {
    let match = null;
    tableCache.forEach(entry => {
      if(!match && deepCompare(params, entry.params)) {
        match = deepCopy(entry);
      }
    });
    return match;
  }

  function loadCache(cache) {
    if(!cache) {
      return;
    }
    rM = cache.rM;
    rdM = cache.rdM;
    rsM = cache.rsM;
    standM = cache.standM;
    hitM = cache.hitM;
    doubleM = cache.doubleM;
    splitM = cache.splitM;
    table = cache.table;
    return cache;
  }

  function cacheTable(params, table) {
    if(!searchTableCache(params)) {
      tableCache.push({
        params: deepCopy(params),
        table: deepCopy(table),
        rM: deepCopy(rM),
        rdM: deepCopy(rdM),
        rsM: deepCopy(rsM),
        standM: deepCopy(standM),
        hitM: deepCopy(hitM),
        doubleM: deepCopy(doubleM),
        splitM: deepCopy(splitM)
      });
    }
  }

  function getTable(player, dealer, table_) {
    if(!table_) {
      table_ = table;
    }
    if(player && dealer) {
      if(dealer === 'A' || dealer === ACE) dealer = DEALER_ACE;
      if(dealer === 10) dealer = DEALER_TEN;
      return table_[TABLE_HANDS.indexOf(player)][DEALER_STATES.indexOf(dealer)];
    } else {
      return table_;
    }
  }

  // P[state][card] = new state
  function progressionMatrix() {
    let a = [];

    HAND_STATES.forEach((hand, i) => {
      a[i] = [];
      let [value, soft] = getHandDetails(hand);
      CARD_STATES.forEach((card, j) => {
        let cardValue = card === ACE ? 1 : card;
        if((hand === DEALER_TEN && card === ACE) || (hand === DEALER_ACE && card === 10)) {
          // Blackjack
          a[i][j] = BLACKJACK;
        } else if(value + cardValue === 21 || (card === ACE && value === 10)) {
          // 21, but not blackjack
          a[i][j] = 21;
        } else if(card === ACE && value < 11) {
          // <11 to soft hand
          a[i][j] = value + ACE;
        } else if(soft && value + cardValue > 21) {
          // Ace goes from 10 -> 1
          a[i][j] = value + cardValue - 10;
        } else if(value + cardValue > 21) {
          // Bust
          a[i][j] = BUST;
        } else {
          // hard -> hard or soft -> soft
          a[i][j] = (value + cardValue) | soft;
        }
      });
    });

    return a;
  }

  // Create transition matrix for the player on hit
  function transitionMatrix(comp) {
    let a = [];
    HAND_STATES.forEach((hand1, i) => {
      a[i] = [];
      HAND_STATES.forEach((hand2, j) => {
        /* If card c brings hand1 to hand2, increase a[i][j] by comp[c]
         * We need to iterate through each card because some hands can be reached from multiple other hands
         * ie. BUST can be reached from 12, 13, 14, etc. */
        a[i][j] = CARD_STATES.reduce((sum, card, c) =>
          sum + comp[c] * (P[i][c] === hand2), 0
        );
      });
    });

    return a;
  }

  function endMatrix(comp, soft17, cards) {
    let states = [];

    // Determine end hands
    let endHands = [17, 18, 19, 20, 21, 18 | SOFT, 19 | SOFT, 20 | SOFT, BLACKJACK, BUST];
    if(!soft17) {
      endHands.push(17 | SOFT);
    }

    DEALER_STATES.forEach((dealer, c) => {
      let shiftedComp;
      if(params.count.system === 'none') {
        shiftedComp = pullCard(comp, c, cards - 1);
      } else if(params.count.system === 'perfect') {
        shiftedComp = copy(comp);
      } else {
        // Using the dealer card, shift the count back to what it was before the dealer card was shown
        // Then using that comp, set shiftedComp to pullCard(comp, c, cards - 1)
        shiftedComp = pullCard(deckComp({
          system: params.count.system,
          indices: params.count.indices,
          tc: params.count.tc + params.count.indices[c] / params.count.decks
        }), c, cards - 1);
      }

      // Progress the state 13 times (maximum size of dealer hand in blackjack)
      states[c] = squishMatrix(progressMState(
        initStateMatrix(dealer, shiftedComp),
        cards - 1,
        endHands,
        13
      ));
    });

    // If peek, weight end states such that blackjack is impossible
    if(params.peek) {
      states.forEach((state, c) => {
        let bjOdds = state[INDEX_BLACKJACK];
        state[INDEX_BLACKJACK] = 0; // Blackjack is impossible
        if(bjOdds > 0) {
           // Weight every other outcome
           states[c] = state.map(endHand => endHand / (1 - bjOdds));
        }
      });
    }

    return states;
  }

  // Calculate player's return by standing
  function standReturns() {
    let a = [];
    HAND_STATES.forEach((hand, i) => {
      a.push([]);
      let value = getHandDetails(hand)[0]; // Player value
      DEALER_STATES.forEach((dealer, j) => {
        let endOdds = endM[j];
        let r = 0; // Return on bet

        endOdds.forEach((odds, k) => {
          let dHand = HAND_STATES[k];
          let v = getHandDetails(dHand)[0]; // Dealer value
          if(hand === BUST) {
            r -= odds;
          } else {
            r += odds * ((value > v) - (value < v));
          }
        });

        a[i].push(r);
      });
    });
    return a;
  }

  function doubleReturns(comp, standM) {
    let a = [];
    let T = transitionMatrix(comp);
    HAND_STATES.forEach((hand, i) => {
      if(hand < params.double.min) {
        // -Infinity return to represent impossible to double
        a[i] = fillArray(-Infinity, T[i].length);
      } else {
        // Matrix multiply T[i] with standM to get return on stand after hitting exactly one card
        // Double the return because we increased our bet
        a[i] = mmultiply(T[i], transpose(standM)).map(x => 2*x);
      }
    });
    return a;
  }

  // Calculate player's return by hitting
  function hitReturns(comp, i, r) {
    return dot(transitionMatrix(comp)[i], r);
  }

  // Calculate player's return by splitting
  /*
   * i is the index of the hand after split
   * ie. if splitting 8s, HAND_STATES[i] = 8
   * r is a vector with the return of each hand for the current dealer card
   * If double after split is allowed, r will be taken from rdM
   * If double after split is not allowed, r will be taken from rM
   * depth=n means we have split n times to reach this hand
   * */
  function splitReturns(comp, cards, i, r, depth) {
    let hand = HAND_STATES[i];
    if(!depth) {
      depth = 1;
    }

    let c = CARD_STATES.indexOf(hand); // Card Index
    // hitReturns is a vector representing the returns of each possible hand after the split
    // ie. pair of 8s -> [return of 10, return of 11, ..., return of soft 19]
    let hitReturns = vmultiply(transitionMatrix(comp)[i], r);
    /*
     * The log2(maxHands) is an approximation that will always give a return below the actual return
     * It assumes that a hand can be split only to a depth of log2(maxHands),
     * however a hand can be split to a further depth if other hands aren't split
     * This approximation is only deals with multiple resplits,
     * and because this is rare, this approximation doesn't affect the return much
     * */
    if(
      params.split.resplit &&
      depth <= Math.log2(params.split.maxHands) &&
      !(hand === ACE && !params.split.resplitAces)
    ) {
      // Important that we shift the comp, or else ten splitting will appear much too advantageous
      comp = pullCard(comp, c, cards.cards);
      cards.cards--;
      // comp[c] is the chance we are able to resplit again
      if(comp[c] > 0) {
        let k = pairIndex(i);
        hitReturns[k] = comp[c] * splitReturns(comp, cards, i, r, depth + 1);
      }
    }
    // Double the return because we increased our bet
    return 2 * vtotal(hitReturns);
  }

  /*-- Move calculation --*/

  function bestMove(comp, hand, dealer, pair) {
    let i = HAND_STATES.indexOf(hand),
        j = DEALER_STATES.indexOf(dealer);

    // Calculate return by splitting
    let split = -Infinity;
    if(pair) {
      let r;
      if(params.split.oneCardAfterAce && hand === ACE) {
        r = standM;
      } else if(params.split.surrender) {
        r = rsM;
      } else if(params.split.double) {
        r = rdM;
      } else {
        r = rM;
      }
      split = splitReturns(
        copy(comp),
        {cards: params.count.cards},
        i,
        transpose(r)[j]
      );
      splitM[CARD_STATES.indexOf(hand)][j] = split;
      hand = pairState(hand); // Hand before split
      i = HAND_STATES.indexOf(hand);
    }

    let stand = standM[i][j],
        double = doubleM[i][j],
        hit;

    // Calculate return by hitting
    hitM[i][j] = hit = hitReturns(comp, i, transpose(params.double.anytime ? rdM : rM)[j]);

    // Determine return on surrender
    let sur = -Infinity;
    if(params.surrender === 'late' || (params.surrender === 'early' && !params.peek)) {
      sur = -.5;
    } else if(params.surrender === 'early') {
      let bjOdds = getBJOdds(comp, dealer);
      sur = (bjOdds - .5) / (1 - bjOdds);
    }

    // Compare returns of each move
    let best;
    if(split > hit && split > stand && split > double) {
      best = new Best('P', split);
    } else if(hit > stand && hit > double) {
      best = new Best('H', hit);
    } else if(stand > double) {
      best = new Best('S', stand);
    } else if(hit > stand) {
      best = new Best('D', double);
    } else {
      best = new Best('d', double);
    }

    // Determine if surrendering is best
    if(sur > best.expected) {
      best.expected = sur;
      best.surrender = true;
    }

    function Best(action, expected) {
      this.action = action;
      this.expected = expected; // Expected return with surrender
      this.expectedNS = expected; // Expected return with no surrender
      this.surrender = false;
    }

    return best;
  }

  function getBJOdds(comp, dealer) {
    if(dealer === DEALER_ACE) dealer = ACE;
    if(dealer === DEALER_TEN) dealer = 10;
    let tempComp = pullCard(comp, CARD_STATES.indexOf(dealer), params.count.decks * 52);
    if(dealer === ACE) {
      return tempComp[CARD_STATES.indexOf(10)];
    } else if(dealer === 10) {
      return tempComp[CARD_STATES.indexOf(ACE)];
    }
    return 0;
  }

  // Chance that we draw each card
  function deckComp(count) {
    // Composition for a deck with no cards pulled
    let baseComp = (() => {
      let a = fillArray(CARD_ODDS, 10);
      a[8] = TEN_ODDS;
      return a;
    })();
    if(count.system === 'ko') {
      return koComp(count);
    } else if(count.system === 'none') {
      return baseComp;
    } else {
      return countingComp(count);
    }
  }

  function countingComp(count) {
    // A perfect counting system works by keeping track of comp
    if(count.system === 'perfect') {
      return count.comp;
    }

    let comp = [];
    /* Each counting system associates each card with an index
     * Sum up the positive and negatives indices
     * ie. for hilo, positive and negative are both 5
     * */
    let [positives, negatives] = count.indices.reduce((acc, index, i) => {
      let multiplier = CARD_STATES[i] === 10 ? 4 : 1;
      if(index > 0) {
        acc[0] += multiplier * index;
      } else if(index < 0) {
        acc[1] -= multiplier * index;
      }
      return acc;
    }, [0, 0]);

    if(negatives === positives) {
      // Balanced system
      comp = count.indices.map((index, i) => {
        let multiplier = CARD_STATES[i] === 10 ? 4 : 1;
        if(index === 0) {
          return CARD_ODDS;
        } else if(index > 0) {
          return multiplier * (CARD_ODDS - count.tc / 104 / positives);
        } else {
          return multiplier * (CARD_ODDS + count.tc / 104 / positives);
        }
      });
    }
    // TODO: Support unbalanced systems

    return comp;
  }

  function koComp(count) {
    const TEN_HIGH_RATIO = 1 - ACE_HIGH_RATIO;

    let comp = [];
    loop(0, 10, i => {
      if(i === 9) {
        // Ace
        comp.push(1.1 * CARD_ODDS + ACE_HIGH_RATIO * count.tc / 104);
      } else if(i === 8) {
        // 10
        comp.push(1.1 * TEN_ODDS + TEN_HIGH_RATIO * count.tc / 104);
      } else if(i <= 5) {
        // 2-7
        comp.push(5.5/6 * CARD_ODDS - count.tc / 104 / 6);
      } else {
        // 8-9
        comp.push(CARD_ODDS);
      }
    });
    return comp;
  }

  // Shift the deck distribution after pulling one card
  function pullCard(state, c, cards) {
    if(state[c] === 0) return zeroes([state.length]);

    // If no cards left in shoe, pull from a newly shuffled shoe
    if(cards <= 0) {
      return fillArray(i => {
        if(CARD_STATES[i] === 10) {
          return TEN_ODDS;
        } else {
          return CARD_ODDS;
        }
      }, 10);
    }

    let newState = [];
    let total = vtotal(state);
    state = vnormalize(state);
    CARD_STATES.forEach((card, c2) => {
      if(c2 === c) {
        newState[c2] = (cards * state[c2] - 1) / (cards - 1);
      } else {
        newState[c2] = (cards * state[c2]) / (cards - 1);
      }
      if(newState[c2] < 0) {
        newState[c2] = 0;
      }
    });
    newState = vnormalize(newState);
    return vscale(newState, total);
  }

  function createHand(value, soft, pair) {
    if(pair) value |= PAIR;
    if(soft) value |= SOFT;
    return value;
  }

  // Returns [value, soft, pair]
  function getHandDetails(hand) {
    if(hand >= 0) {
      return [hand & VALUE_MASK, hand & SOFT, hand & PAIR];
    } else if(hand === BLACKJACK) {
      return [Infinity, 0, 0];
    } else if(hand === BUST) {
      return [-Infinity, 0, 0];
    } else if(hand === DEALER_TEN) {
      return [10, 0, 0];
    } else if(hand === DEALER_ACE) {
      return [11, SOFT, 0];
    } else {
      return [0, 0, 0];
    }
  }

  function shuffle(comp, cards) {
    cards.cards = Math.round(52 * params.count.decks);
    if(cards.count !== undefined) {
      cards.count = 0;
      cards.trueCount = 0;
    }
    return copy(comp);
  }

  /*-- General Utility Functions --*/

  // mState[i][j]: the chance of player having hand i and getting card j on next hit
  // cards: number of cards
  // endHands[i]: don't hit on hand i
  function progressMState(mState, cards, endHands, reps) {
    let newState = zeroes([HAND_STATES.length, CARD_STATES.length]);

    if(reps) {
      for(let i = 0; i < reps; i++) {
        mState = progressMState(mState, cards--, endHands);
      }
      return mState;
    }

    HAND_STATES.forEach((postHit, I) => {
      HAND_STATES.forEach((preHit, i) => {
        // preHit: hand before hitting
        // postHit: hand after hitting

        if(endHands.includes(postHit) && i === I) {
          // If postHit is an end hand, we don't draw any cards
          return newState[I] = vsum(newState[I], mState[I]);
        }

        let [value, soft] = getHandDetails(preHit);
        CARD_STATES.forEach((card, j) => {
          let preHitChance = vtotal(mState[i]); // Chance that player has hand preHit
          // If card causes preHit to transition to postHit (eg. card=2, preHit=8, postHit=10)
          if(
            preHitChance > 0 &&
            !endHands.includes(preHit) &&
            P[i][j] === postHit
          ) {
            // Shift the composition as if we pulled card j
            let newComp = pullCard(mState[i], j, cards);
            // Add that deck distribution to (state I)
            newState[I] = vsum(newState[I], vscale(newComp, mState[i][j] / preHitChance));
          }
        });
      });
    });

    return newState;
  }

  function deepCompare(obj1, obj2) {
    if(
      obj1 !== undefined && obj2 === undefined ||
      obj1 === undefined && obj2 !== undefined
    ) return false;

    let propertiesMatch = true;
    Object.keys(obj1).forEach((key) => {
      if(typeof obj1[key] === 'object') {
        if(!obj2.hasOwnProperty(key) || !deepCompare(obj1[key], obj2[key])) {
          propertiesMatch = false;
        }
      } else if(obj1[key] !== obj2[key]) {
        propertiesMatch = false;
      }
    });

    // Verify that keys match
    let keysMatch = true;
    Object.keys(obj2).forEach((key) => {
      if(!obj1.hasOwnProperty(key)) {
        keysMatch = false;
      }
    });

    return propertiesMatch && keysMatch;
  }

  function pairState(state) {
    return state = state === ACE ? 44 : state * 2;
  }

  function pairIndex(i) {
    return HAND_STATES.indexOf(pairState(HAND_STATES[i]));
  }

  /*-- Logic Utility Functions --*/

  function loop(start, end, f) {
    for(let i = start; i < end; i++) {
      f(i);
    }
  }

  // Like loop but in reverse
  function rLoop(start, end, f) {
    for(let i = start; i > end; i--) {
      f(i);
    }
  }

  function fillArray(value, len) {
    let a = [];
    if(typeof value === 'function') {
      loop(0, len, () => a.push(value()));
    } else {
      loop(0, len, () => a.push(value));
    }
    return a;
  }

  function zeroes(dims) {
    if(dims.length === 1) {
      return fillArray(0, dims[0]);
    } else {
      return fillArray(zeroes.bind(null, dims.slice(1)), dims[0]);
    }
  }

  function copy(v) {
    let u = [];
    for(let i = 0; i < v.length; i++) {
      u[i] = v[i];
    }
    return u;
  }

  function deepCopy(obj) {
    if(typeof obj !== 'object') return obj;
    if(obj.length !== undefined) {
      let copy = [];
      obj.forEach(elem => {
        copy.push(deepCopy(elem));
      });
      return copy;
    } else {
      let copy = {};
      Object.keys(obj).forEach(key => {
        copy[key] = deepCopy(obj[key]);
      });
      return copy;
    }
  }

  /*-- Mathy Utility Functions --*/

  function vnormalize(u) {
    return vscale(u, 1/vtotal(u));
  }

  function vtotal(u) {
    return u.reduce((r, x) => r + x, 0);
  }

  function vsum(u, v) {
    if(typeof v === 'object') {
      return u.map((x, i) => x + v[i]);
    } else {
      return u.map((x, i) => x + v);
    }
  }

  function vscale(u, c) {
    return u.map((x, i) => x * c);
  }

  function dot(u, v) {
    return u.reduce((r, x, i) => r + x * v[i], 0);
  }

  function vmultiply(u, v) {
    return u.map((x, i) => x * v[i]);
  }

  function transpose(a) {
    let b = [];
    loop(0, a[0].length, j => {
      b.push([]);
      loop(0, a.length, i => b[j].push(a[i][j]));
    });
    return b;
  }

  // mmultiply(vec, matrix) or mmultiply(matrix, matrix)
  function mmultiply(v, a) {
    if(v[0] && v[0].length !== undefined) {
      let c = [];
      a = transpose(a);
      // v times a
      v.forEach((row, i) => {
        c.push([]);
        a.forEach((col, j) => c[i].push(dot(row, col)));
      });
      return c;
    } else {
      let w = [];
      for(let u of a) {
        w.push(dot(v, u));
      }
      return w;
    }
  }

  // Raise a matrix to a positive integer power (n>0)
  function mpower(a, n) {
    let p = a;
    loop(1, n, () => p = mmultiply(p, a));
    return p;
  }

  // Picks random index from weighted probability array
  function pickFromArray(a) {
    let x = Math.random() * vtotal(a),
        i = -1;
    while(x > 0) {
      x -= a[++i];
    }
    return i;
  }

  // 2d matrtix -> 1d vector
  function squishMatrix(m) {
    let v = [];
    HAND_STATES.forEach((hand, i) => {
      v[i] = 0;
      CARD_STATES.forEach((card, j) => {
        v[i] += m[i][j];
      });
    });
    return v;
  }

  function initStateMatrix(hand, comp) {
    let stateMatrix = zeroes([HAND_STATES.length, CARD_STATES.length]);
    stateMatrix[HAND_STATES.indexOf(hand)] = copy(comp);
    return stateMatrix;
  }

  callback(); // Callback
}

let jackfish, sim;
self.addEventListener('message', e => {
 let f = e.data[0];
 let args = e.data[1];
 let callback = messageCallback.bind(null, f);
 switch(f) {
   case 'Constructor':
    jackfish = new Jackfish(args[0], callback);
    break;
  case 'createSimulation':
    sim = jackfish.createSimulation(args[0]);
    break;
  case 'updateSimulation':
    sim.config(args[0]);
    break;
  case 'runSimulation':
    sim.run(callback);
    callback(sim.get());
    break;
  case 'clearSimulation':
    sim.clear();
    callback(sim.get());
    break;
  case 'getSimulation':
    callback(sim.get());
    break;
  case 'stopSimulation':
    sim.stop();
    callback(sim.get());
    break;
  case 'setParams':
    if(jackfish === undefined) {
      jackfish = new Jackfish(args[0], callback);
    }
    jackfish.setParams(callback, args[0]);
    if(sim) {
      sim.config(sim.getOptions());
    }
    break;
  case 'getParams':
    jackfish.getParams(callback);
    break;
  case 'getTable':
    jackfish.getTable(callback);
    break;
  case 'doAll':
    jackfish.doAll(callback);
    break;
  case 'makeMatrices':
    jackfish.makeMatrices(callback);
    break;
  case 'makeTable':
    jackfish.makeTable(callback);
    break;
  case 'takeInsurance':
    jackfish.takeInsurance(callback);
    break;
  case 'getEdge':
    jackfish.getEdge(callback, args[0]);
    break;
 }
});

function messageCallback(f, data) {
  postMessage([f, data]);
}
