const path = require('path');
require('dotenv').config({ path: path.join(__dirname, '.env') });
const unifiedAgent = require('./src/unified-agent');
const e2bManager = require('./src/e2b-manager');
const fs = require('fs');

// Q39 - Jupiter moons resonance (with ALL options)
const question39 = {
  id: "mmlu-39",
  category: "physics",
  question: `What is the significance of the 1:2:4 resonance in the Jupiter's moons system?

Options:
(A) It causes the three moons to always be in a straight line.
(B) It causes Io to rotate on its axis faster than it orbits Jupiter
(C) It creates a gap with no asteroids between the orbits.
(D) It results in the formation of a visible ring around Jupiter.
(E) It prevents the moons from colliding with each other.
(F) It causes a gravitational pull that affects Jupiter's orbit around the sun.
(G) It results in a uniform gravitational pull on all three moons.
(H) It prevents formation of the ring material into other moons.
(I) It makes the orbit of Io slightly elliptical.`,
  correct_answer: "I"
};

// Q40 - Pluto diameter measurement (with ALL options)
const question40 = {
  id: "mmlu-40",
  category: "physics",
  question: `We were first able to accurately measure the diameter of Pluto from:

Options:
(A) Lunar-based observations made during NASA's Apollo missions
(B) The Voyager 2 spacecraft flyby in the 1980s
(C) Observations made by the Chandra X-ray Observatory
(D) Brightness measurements made during mutual eclipses of Pluto and Charon
(E) The Mars Rover's telescope observations of Pluto
(F) Radar observations made by the Arecibo telescope
(G) Observations made by the Spitzer Space Telescope
(H) Hubble Space Telescope images that resolved Pluto's disk
(I) A New Horizons flyby in the 1990s
(J) Images captured by the Kepler Space Telescope`,
  correct_answer: "D"
};

async function testQuestions() {
  console.log('=== RETESTING Q39 & Q40 WITH COMPLETE OPTIONS ===\n');
  
  // Test Q39
  console.log('Testing Q39 - Jupiter moons resonance');
  console.log('Correct answer: I (It makes the orbit of Io slightly elliptical)\n');
  
  let startTime = Date.now();
  try {
    const result39 = await unifiedAgent.process(question39.question, {}, {
      requestId: 'mmlu-q39-retest-' + Date.now()
    });
    
    const filename39 = `mmlu-q39-retest-${Date.now()}.txt`;
    fs.writeFileSync(filename39, result39.content);
    
    console.log('Q39 Response time:', ((Date.now() - startTime) / 1000).toFixed(1), 's');
    console.log('Response saved to:', filename39);
    console.log('\nLast 600 chars:');
    console.log(result39.content.substring(Math.max(0, result39.content.length - 600)));
    
  } catch (error) {
    console.error('Q39 Error:', error.message);
  }
  
  console.log('\n' + '='.repeat(60) + '\n');
  
  // Test Q40
  console.log('Testing Q40 - Pluto diameter measurement');
  console.log('Correct answer: D (Brightness measurements during mutual eclipses)\n');
  
  startTime = Date.now();
  try {
    const result40 = await unifiedAgent.process(question40.question, {}, {
      requestId: 'mmlu-q40-retest-' + Date.now()
    });
    
    const filename40 = `mmlu-q40-retest-${Date.now()}.txt`;
    fs.writeFileSync(filename40, result40.content);
    
    console.log('Q40 Response time:', ((Date.now() - startTime) / 1000).toFixed(1), 's');
    console.log('Response saved to:', filename40);
    console.log('\nLast 600 chars:');
    console.log(result40.content.substring(Math.max(0, result40.content.length - 600)));
    
  } catch (error) {
    console.error('Q40 Error:', error.message);
  }
  
  await e2bManager.shutdown();
}

testQuestions().catch(console.error);