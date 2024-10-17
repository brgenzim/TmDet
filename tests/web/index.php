<?php

ob_start();

if (($pdbCode = get_cfg_var('pdbCode')) === false) {
    error('Empty PDB CODE');
}
$molstarId = 'molstar_' . $pdbCode;


?><!DOCTYPE html>
<html lang="en">
<head>
    <title>TmDet :: <?= $pdbCode ?></title>
    <link href="/css/molstar.css" rel="stylesheet"/>
    <link href="/css/molstar_viewer.css" rel="stylesheet"/>
    <link href="/wrench-adjustable.svg" rel="icon" type="image/svg+xml" />
</head>
<body>

<div style="width: 50%;">
<div class="title"><?= $pdbCode ?></div>
<div class="molstar-viewer" id="<?= $molstarId ?>"></div>
</div>

<script type="text/javascript" src="/js/tm_molstar.js"></script>

<script type="text/javascript">

    var viewer = new tm_molstar.Viewer("<?= $molstarId ?>", {
        layoutShowControls: false,
        layoutIsExpanded: false,
        viewportShowExpand: true,
        collapseLeftPanel: true
    });

    tm_molstar.loadWithUNITMPMembraneRepresentation(viewer.plugin, {
        structureUrl: "/cif.php",
        regionDescriptorUrl: "/json.php"
    });

</script>
</body>
</html>

<?php
function error(string $message) {
    ob_end_clean();
    header(header: 'HTTP/1.1 500 Internal Server Error', response_code: 500);
    header('Content-type: application/json');
    echo json_encode([ 'error' => $message ], JSON_UNESCAPED_SLASHES);
    exit();
}
?>

<?php ob_end_flush(); ?>
