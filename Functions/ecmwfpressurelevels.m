

% returns the ROUGH ecmwf pressure/altitude levels, depending on the number of
% levels you enter. Default is 137 levels.


function [pres,alt] = ecmwfpressurelevels(IN)

if nargin == 0
    numlevels = 137;
else
    numlevels = IN;
end

switch nargout
    case {0 1}
        pres = getpressurelevels(numlevels);
    case 2
        pres = getpressurelevels(numlevels);
        alt  = getaltitudelevels(numlevels);
end

end

% MODEL PRESSURE LEVEL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pres = getpressurelevels(numlevels)

switch numlevels
    
    case 137
        
        pres = [0.02
            0.031
            0.0467
            0.0683
            0.0975
            0.1361
            0.1861
            0.2499
            0.3299
            0.4288
            0.5496
            0.6952
            0.869
            1.0742
            1.3143
            1.5928
            1.9134
            2.2797
            2.6954
            3.1642
            3.6898
            4.2759
            4.9262
            5.6441
            6.4334
            7.2974
            8.2397
            9.2634
            10.372
            11.5685
            12.8561
            14.2377
            15.7162
            17.2945
            18.9752
            20.761
            22.6543
            24.6577
            26.7735
            29.0039
            31.3512
            33.8174
            36.4047
            39.1149
            41.9493
            44.9082
            47.9915
            51.199
            54.5299
            57.9834
            61.5607
            65.2695
            69.1187
            73.1187
            77.281
            81.6182
            86.145
            90.8774
            95.828
            101.0047
            106.4153
            112.0681
            117.9714
            124.1337
            130.5637
            137.2703
            144.2624
            151.5493
            159.1403
            167.045
            175.2731
            183.8344
            192.7389
            201.9969
            211.6186
            221.6146
            231.9954
            242.7719
            253.9549
            265.5556
            277.5852
            290.0548
            302.9762
            316.3607
            330.2202
            344.5663
            359.4111
            374.7666
            390.645
            407.0583
            424.019
            441.5395
            459.6321
            478.3096
            497.5845
            517.4198
            537.7195
            558.343
            579.1926
            600.1668
            621.1624
            642.0764
            662.8084
            683.262
            703.3467
            722.9795
            742.0855
            760.5996
            778.4661
            795.6396
            812.0847
            827.7756
            842.6959
            856.8376
            870.2004
            882.791
            894.6222
            905.7116
            916.0815
            925.7571
            934.7666
            943.1399
            950.9082
            958.1037
            964.7584
            970.9046
            976.5737
            981.7968
            986.6036
            991.023
            995.0824
            998.8081
            1002.225
            1005.3562
            1008.2239
            1010.8487
            1013.25];
        
    case 91
        
        pres = [0.02
            0.0398
            0.0739
            0.1291
            0.2141
            0.3395
            0.5175
            0.7617
            1.0872
            1.5099
            2.0464
            2.7136
            3.5282
            4.5069
            5.6652
            7.0181
            8.5795
            10.3617
            12.3759
            14.6316
            17.1371
            19.8987
            22.9216
            26.209
            29.763
            33.5843
            37.672
            42.0242
            46.6378
            51.5086
            56.6316
            61.9984
            67.5973
            73.415
            79.4434
            85.7016
            92.2162
            99.0182
            106.1445
            113.6382
            121.5502
            129.9403
            138.8558
            148.326
            158.3816
            169.0545
            180.3786
            192.3889
            205.1222
            218.6172
            232.914
            248.0547
            264.0833
            281.0456
            298.9895
            317.9651
            338.0245
            359.2221
            381.6144
            405.2606
            430.2069
            456.4813
            483.8505
            512.0662
            540.8577
            569.9401
            599.031
            627.9668
            656.6129
            684.8491
            712.5573
            739.5739
            765.7697
            791.0376
            815.2774
            838.3507
            860.1516
            880.608
            899.6602
            917.2205
            933.2247
            947.6584
            960.5245
            971.8169
            981.5301
            989.7322
            996.8732
            1002.8013
            1007.4431
            1010.8487
            1013.25];
        
    case 60
        pres = [0.2
            0.3843
            0.6365
            0.9564
            1.3448
            1.8058
            2.3478
            2.985
            3.7397
            4.6462
            5.7565
            7.1322
            8.8366
            10.9483
            13.5647
            16.8064
            20.8227
            25.7989
            31.9642
            39.6029
            49.0671
            60.1802
            73.0663
            87.7274
            104.2288
            122.6137
            142.9017
            165.0886
            189.1466
            215.0251
            242.6523
            272.0593
            303.2174
            336.0439
            370.4072
            406.1328
            443.0086
            480.7907
            519.2093
            557.9734
            596.7774
            635.306
            673.2403
            710.2627
            746.0635
            780.3455
            812.8303
            843.2634
            871.4203
            897.1118
            920.1893
            940.5511
            958.1477
            972.9868
            985.1399
            994.7472
            1002.0236
            1007.2639
            1010.8487
            1013.25];
        
end
end

% MODEL ALTITUDE LEVEL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alt = getaltitudelevels(numlevels)

switch numlevels
    
    case 137
        
        alt = [80301.65
            74584.91
            71918.79
            69365.77
            66906.53
            64537.43
            62254.39
            60053.46
            57930.78
            55882.68
            53905.62
            51996.21
            50159.36
            48413.94
            46756.98
            45199.69
            43738.55
            42364.93
            41071.2
            39850.56
            38696.94
            37604.95
            36569.72
            35586.89
            34652.52
            33763.05
            32915.27
            32106.57
            31330.96
            30584.71
            29866.09
            29173.5
            28505.47
            27860.64
            27237.73
            26635.56
            26053.04
            25489.15
            24942.93
            24413.5
            23900.02
            23401.71
            22917.85
            22447.75
            21990.82
            21546.62
            21114.77
            20694.9
            20286.66
            19889.88
            19503.09
            19125.61
            18756.34
            18394.25
            18038.35
            17687.77
            17341.62
            16999.08
            16659.55
            16322.83
            15988.88
            15657.7
            15329.24
            15003.5
            14680.44
            14360.05
            14042.3
            13727.18
            13414.65
            13104.7
            12797.3
            12492.44
            12190.1
            11890.24
            11592.86
            11297.93
            11005.69
            10714.22
            10422.64
            10130.98
            9839.26
            9547.49
            9255.7
            8963.9
            8672.11
            8380.36
            8088.67
            7797.04
            7505.51
            7214.09
            6922.8
            6631.66
            6340.68
            6049.89
            5759.3
            5469.3
            5180.98
            4896.02
            4615.92
            4341.73
            4074.41
            3814.82
            3563.69
            3321.67
            3089.25
            2866.83
            2654.69
            2452.99
            2261.8
            2081.09
            1910.76
            1750.63
            1600.44
            1459.91
            1328.7
            1206.44
            1092.73
            987.15
            889.29
            798.72
            715.02
            637.76
            566.54
            500.95
            440.61
            385.16
            334.24
            287.52
            244.69
            205.44
            169.51
            136.62
            106.54
            79.04
            53.92
            30.96
            10];
        
    case 91
        
        alt = [79302.74
            72744.48
            68689.83
            64847.81
            61203.87
            57748.06
            54470.24
            51360.33
            48442.75
            45759.1
            43332.12
            41136.44
            39141.66
            37322.4
            35657.28
            34128.22
            32719.8
            31417.81
            30201.49
            29061.56
            27991.33
            26984.85
            26036.84
            25142.55
            24297.75
            23498.64
            22741.78
            22024.08
            21342.73
            20695.19
            20079.13
            19492.24
            18931.56
            18396.14
            17884.58
            17394.28
            16921.83
            16464.09
            16018.24
            15581.7
            15152.05
            14727.04
            14305
            13885.42
            13468.26
            13053.52
            12641.18
            12231.23
            11823.65
            11418.44
            11015.57
            10613.07
            10207.82
            9800.09
            9389.89
            8977.23
            8562.14
            8144.61
            7724.68
            7302.35
            6877.76
            6451.33
            6025.61
            5604.74
            5192.42
            4791.93
            4406.04
            4036.32
            3683.17
            3346.69
            3026.84
            2723.84
            2437.96
            2169.15
            1917.26
            1682.35
            1464.53
            1263.63
            1079.33
            911.48
            759.99
            624.53
            504.53
            399.5
            309.03
            232.49
            167.39
            112.26
            67.88
            34.22
            10];
        
    case 60
        
        alt = [64947.31
            57350.89
            53125.16
            49623.45
            46709.58
            44259.61
            42156.09
            40294.73
            38600.97
            37018.28
            35500.64
            34017.99
            32561.15
            31127.02
            29702.74
            28287.36
            26880.84
            25483.11
            24094.12
            22713.82
            21342.15
            20014.52
            18755.37
            17563.62
            16440.2
            15381.19
            14382.89
            13441.79
            12554.62
            11718.27
            10930.02
            10175.26
            9444.61
            8737.51
            8054.2
            7395.38
            6761.89
            6154.67
            5574.58
            5022.43
            4498.91
            4004.59
            3539.96
            3105.35
            2701
            2327.04
            1983.49
            1670.26
            1387.12
            1133.73
            909.57
            713.97
            546.06
            404.72
            288.55
            195.85
            124.48
            71.89
            34.97
            10];
        
        
end

end








