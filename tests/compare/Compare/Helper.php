<?php

namespace Tmdet\Compare;

class Helper {

    public static function clearSeq(string $raw) : string {
        $ret = "";
        for($i=0; $i<strlen($raw); $i++) {
            if (strchr(" \t\n",$raw[$i]) === false) {
                $ret .= ($raw[$i]=='?'?'X':$raw[$i]);
            }
        }
        return $ret;
    }
}