<?php

$differences = json_decode(json: file_get_contents($argv[1]), associative: true);
$entry = $differences['detailsByCodes'][$argv[2]];
echo json_encode($entry, JSON_PRETTY_PRINT) . PHP_EOL;