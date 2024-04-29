<?php

namespace Unitmp\TmdetTest;

use PHPUnit\Framework\Attributes\Group;
use PHPUnit\Framework\TestCase;
use Unitmp\TmdetTest\Services\DsspRunner;
use Unitmp\TmdetTest\Services\TmDetDsspRunner;

class TmDetDsspTest extends TestCase {

    #[Group('dssp')]
    public function test_1afo_dssp() {
        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner('1AFO');
        $dsspRunner = DsspRunner::createRunner('1AFO');
        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        $dsspSuccess = $dsspRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
    }

    #[Group('dssp')]
    public function test_7rh7_dssp() {
        // Arrange
        $tmdetRunner = TmDetDsspRunner::createRunner('7RH7');
        $dsspRunner = DsspRunner::createRunner('7RH7');
        // Act
        $tmdetSuccess = $tmdetRunner->exec();
        $dsspSuccess = $dsspRunner->exec();
        // Assert
        $this->assertTrue($tmdetSuccess);
        $this->assertTrue($dsspSuccess);
        $this->assertEquals($dsspRunner->dssps, $tmdetRunner->dssps);
    }

}
