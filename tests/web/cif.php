<?php

ob_start();

$cifFile = getPath();
$fileName = basename($cifFile);
$data = file_get_contents($cifFile);

$isGzipped = str_ends_with($cifFile, '.gz');
$isZipped = str_ends_with($cifFile, '.zip');
if ($isGzipped || $isZipped) {
    header('Content-Type: application/zip');
    if ($isGzipped) {
        header('Content-Encoding: gzip');
    }
    if ($isZipped) {
        header('Content-Encoding: deflate');
    }
} else {
    header('Content-Type: text/plain');
}

header("Content-Disposition: attachment; filename=\"$fileName\"");
header('Content-Length: ' . strlen($data));

echo $data;
ob_end_flush();

//
// Utils
//

function getPath(): string {

    if (empty($cifPath = get_cfg_var('cifFile'))) {
        error('Empty cifFile variable');
    }

    if (!file_exists($cifPath)) {
        error("File not found: $cifPath", 404);
    }

    return $cifPath;
}

function error(string $message, int $code = 500) {
    ob_end_clean();
    if ($code == 500) {
        header(header: 'HTTP/1.1 500 Internal Server Error', response_code: $code);
    } else {
        header(header: "HTTP/1.1 $code Error", response_code: $code);
    }
    header('Content-Type: application/json');
    echo json_encode([ 'error' => $message ], JSON_UNESCAPED_SLASHES);
    exit();
}

