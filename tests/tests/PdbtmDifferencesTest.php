<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

class PdbtmDifferencesTest extends TestCase {

    const PDB_ENT_DIR = FileSystem::PDB_ENT_ZFS_DIR;
    const MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT = 2;

    #[Group('tmdet')]
    #[DataProvider('getEntPathsOfTmps')]
    public function test_pdbtm_manual_updates(string $entPath) {
        // Arrange
        $tmdetRunner = Tmdet30Runner::createRunner($entPath);
        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertEquals($tmdetRunner->oldData, $tmdetRunner->newData);
    }

    //
    // Util functions
    //


    public static function getEntPathsOfTmps(): array {
        $paths = [];
        foreach (explode("\n", trim(file_get_contents('../data/tm-codes.txt'))) as $code) {
            $fileName = "$code.pdb.gz";
            $subDir = '/' . $code[1] . $code[2] . '/';
            $paths[$fileName] = [ static::PDB_ENT_DIR . $subDir . $fileName  ];
        }
        return $paths;
        //return array_slice($paths, 0, 10);
    }

}
