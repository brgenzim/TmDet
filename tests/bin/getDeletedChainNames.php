<?php

$differences = json_decode(json: file_get_contents($argv[1]), associative: true);

$chainNames = [];
$detailsByCode = $differences['detailsByCodes'];
// $detailsByCode = array_slice($detailsByCode, 0, 100);
foreach ($detailsByCode as $entry) {
    $code = $entry['code'];
    $details = $entry['details'];
    if (!array_key_exists('oldData', $details) || !array_key_exists('deletedChains', $details['oldData'])) {
        continue;
    }
    foreach ($details['oldData']['deletedChains'] as $chain) {
        $polymers = $entry['polymerChains'];
        if (!array_key_exists($chain, $polymers)) {
            continue;
        }
        if (array_key_exists($chain, $details['newData']['deletedChains'])) {
            continue;
        }

        $name = $polymers[$chain];
        if (empty(trim($name))) {
            // var_dump([ $code, $entry['polymerChains'] ]);
        }
        if (!array_key_exists($name, $chainNames)) {
            $chainNames[$name] = [ $code ];
        } elseif (array_search($code, $chainNames[$name]) === false) {
            array_push($chainNames[$name], $code);
        }
    }
}

ksort($chainNames);

echo json_encode($chainNames, JSON_PRETTY_PRINT|JSON_UNESCAPED_SLASHES) . PHP_EOL;
