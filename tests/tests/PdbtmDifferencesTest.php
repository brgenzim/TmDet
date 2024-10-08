<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\DataProvider;
use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Constants\FileSystem;
use Unitmp\TmdetTest\Services\PdbtmComparator;
use Unitmp\TmdetTest\Services\Tmdet30Runner;

class PdbtmDifferencesTest extends TestCase {

    const PDB_ENT_DIR = FileSystem::PDB_ENT_ZFS_DIR;
    const PDBTM_ALL_XML_FILE = FileSystem::PDBTM_ALL_XML_FILE;
    const MAX_ALLOWED_LEVENSHTEIN_DISTANCE_IN_PERCENT = 40;

    #[Group('tmdet')]
    #[DataProvider('getEntPathsOfTmps')]
    public function test_pdbtm_manual_updates(string $entPath) {
        // Arrange
        $tmdetRunner = Tmdet30Runner::createRunner($entPath);
        if ($tmdetRunner->oldData['tmType'] !== 'Tm_Alpha') {
            printf("Skipping test (not Tm_Alpha): %s\n", $tmdetRunner->oldTmdetFile);
            $this->markTestSkipped('Not Tm_Alpha ' . $entPath);
        }

        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $result = PdbtmComparator::compareTmdetData($tmdetRunner);
        // if there is any difference, the first one will be the message
        $messages = $result['messages'];
        $this->assertEmpty($messages, empty($messages) ? 'no error' : $messages[0]);
    }

    //
    // Util functions
    //

    public static function getEntPathsOfTmps(): array {
        $codes = Tmdet30Runner::getPdbtmCodes(static::PDBTM_ALL_XML_FILE);

        $paths = [];
        foreach ($codes as $code => $item) {
            $subDir = '/' . $code[1] . $code[2] . '/';
            $fileName = "pdb$code.ent.gz";
            $fullPath = static::PDB_ENT_DIR . $subDir . $fileName;
            $paths[ $fileName ] = [ $fullPath ];
        }
        ksort($paths);
        //return array_slice($paths, 0, 100);
        return $paths;
    }

    // public static function getEntPathsOfTmps(): array {
    //     $codes = explode(PHP_EOL, trim(file_get_contents('issue_943/tmdet_xml_differences_differences_v4.2.failures.different.chains.txt')));
    //     $codes = array_slice($codes, 0, 100);

    //     $paths = [];
    //     foreach ($codes as $code) {
    //         $subDir = '/' . $code[1] . $code[2] . '/';
    //         $fileName = "pdb$code.ent.gz";
    //         $fullPath = static::PDB_ENT_DIR . $subDir . $fileName;
    //         $paths[ $fileName ] = [ $fullPath ];
    //     }
    //     ksort($paths);
    //     return $paths;
    // }

}
