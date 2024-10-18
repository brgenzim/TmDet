<?php

ob_start();

$radius = get_cfg_var('radius');
$thickness = get_cfg_var('thickness');
if (empty($radius) || empty($thickness)) {
    error("'radius' and 'thickness' are mandatory config variables");
}
if (empty($pdbCode = get_cfg_var('pdbCode'))) {
    error("Empty 'pdbCode' config variable");
}
$radius = floatval($radius);
$thickness = floatval($thickness);

// generate response body
$data = getJson(pdbCode: $pdbCode, radius: 14, thickness: 30);

header('Content-Type: application/json');
header('Content-Length: ' . strlen($data));

echo $data;


function error(string $message, int $code = 500) {

    ob_end_clean();
    if ($code == 500) {
        header(header: 'HTTP/1.1 500 Internal Server Error', response_code: $code);
        header('Content-type: application/json');
    } else {
        header(header: "HTTP/1.1 $code Error", response_code: $code);
    }
    echo json_encode([ 'error' => $message ], JSON_UNESCAPED_SLASHES);
    exit();
}

function getJson(string $pdbCode, float $radius, float $thickness): string {

    $annotation = [
        "data_resource" => "TmDet Debug Web App",
        "pdb_id" => $pdbCode,
        "includes_het_groups" => false,
        "chains" => [],
        "sites" => [],
        "additional_entry_annotations" => [
            "membrane" => [
                "normal" => [ "x" => 0, "y" => 0, "z" => $thickness / 2.0 ],
                "radius" => $radius,
                "transformation_matrix" => [
                    "rowx" => [ "x" => 1, "y" => 0, "z" => 0, "t" => 0 ],
                    "rowy" => [ "x" => 0, "y" => 1, "z" => 0, "t" => 0 ],
                    "rowz" => [ "x" => 0, "y" => 0, "z" => 1, "t" => 0 ]
                ]
            ]
        ]
    ];
    return json_encode($annotation);
}
