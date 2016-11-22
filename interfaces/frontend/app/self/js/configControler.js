/*
* @Author: jaume.bonet
* @Date:   2016-11-03 14:23:09
* @Last Modified by:   jaume.bonet
* @Last Modified time: 2016-11-11 16:47:38
*/

'use strict';

// var fixMaxLink = function( reference ){
//     return parseFloat((Math.sqrt(2 * Math.pow(reference, 2)) + 2.0).toFixed(2));
// }
// var fillConfig = function( ) {
//     var c = { name: null,
//               hLength: 13, eLength: 7,
//               heDist: 11, hhDist: 12, epDist: 5, etDist: 5,
//               maxlink: 0, keeplinkratio: false,
//               interface: {
//                 host: 'localhost', port: 8080,
//                 lockedCanvas: false, proportion: 10,
//                 helix:        {strokeColor: '#000000', strokeWidth: 1, fillColor: '#ccccff'},
//                 beta:         {strokeColor: '#000000', strokeWidth: 1, fillColor: '#ffcccc'},
//                 cross:        {strokeColor: '#000000', strokeWidth: 1, fillColor: '#ccffcc'},
//                 placeholders: {strokeColor: '#e6b0aa', strokeWidth: 1, fillColor: '#fdfefe'},
//                 motif:        {strokeColor: '#000000', strokeWidth: 1, fillColor: '#ffffff'}
//              }
//         };
//     c.maxlink = fixMaxLink(c.hhDist);
//     return c;
// }
angular.module('topobuilder')
.controller('ConfigCtrl', function($scope, $http, CONFIG, SSE, MOTIFS, configFUNC) {
    var self = this;

    $scope.config  = CONFIG;
    // $scope.config  = fillConfig();
    $scope.motifs  = MOTIFS;
    $scope.config.maxlink = parseFloat((Math.sqrt(2 * Math.pow($scope.config.hhDist, 2)) + 2.0).toFixed(2));
    self.default = _.clone($scope.config);
    $scope.info    = true;
    self.infile = null;

    $scope.toggleInfo = function() {
        $scope.info = !$scope.info;
    };
    $scope.canvas2Layers = function() {
        return _.map(SSE.list, function(p){return p.tpb});
    }

    $scope.project2JSON = function() {
        var project = {config: _.clone($scope.config),
                       layers: $scope.canvas2Layers(),
                       motifs: _.clone($scope.motifs.desc)};
        return project;
    }
    self.JSON2project = function( project ) {
        Object.assign($scope.config, project.config)
        //TODO the rest
    }

    self.loadJob = function() {
        console.log(self.infile);
    };

    self.saveJob = function() {
        var formData = configFUNC.postForm('save')
        formData.append('content', JSON.stringify($scope.project2JSON()));

        $http({
            method: 'POST', url: configFUNC.connection() + '/query',
            headers: { 'Content-Type': undefined }, data: formData,
            transformRequest: angular.identity
        })
        .success(function (data) {
            console.log(data)
            if (data.success == 'ok') { swal("Saved!", "", "success"); }
            else if (data.success == 'ko') { swal("Error!", data.msg, "error"); }
            else { swal("Error!", "Unknown error has occurred.", "error"); }

        })
        .error(function (data, status) {
            console.log(data); console.log(status);
            swal("Error!", "http error: " + status, "error");
        });
    };

    self.fetchJob = function() {
        
    };

    //Link Correlation
    $scope.updateLink = function() {
        if ($scope.config.keeplinkratio) {
            $scope.config.maxlink = parseFloat((Math.sqrt(2 * Math.pow($scope.config.hhDist, 2)) + 2.0).toFixed(2));
        }
    }

    // Tab Management & Visualization
    self.tab = 0;
    self.setTab = function( value ) {
        if (!self.noProjectID()) self.tab = value;
    };
    self.evalTab = function( value ) {
        if (self.noProjectID()) return false;
        return value == self.tab;
    };
    self.noProjectID = function() { return !configFUNC.projectHasStart(); };

    // Host & Ports
    self.isLocalhost = function() {
        return _.contains( ['localhost', '127.0.0.1'], $scope.config.interface.host )
    }
    self.testConnection = function() {
        $http({
            method: 'GET', url: configFUNC.connection(),
        })
        .success(function (data) {
            if (data.services && _.contains(Object.keys(data.services), 'topobuilder')) {
                swal("Connection!", "The selected host/port has the topobuilder service.", "success");
            }
            else {
                swal("Connection Error!", "But the selected host/port does not have the topobuilder service.", "error");
            }
        })
        .error(function (data, status) {
            swal("Connection Error!", "The host/port returns an error.", "error");
        });
    }

    //Default Storage Management
    self.name2default = function() {
        self.default.name = _.clone($scope.config.name);
        self.default.pin  = _.clone($scope.config.pin);
        self.default.date = _.clone($scope.config.date);
    };
    self.resetConfig = function() {
        var pairs = _.pairs(self.default)
        for (var i = pairs.length - 1; i >= 0; i--) {
            $scope.config[pairs[i][0]] = pairs[i][1];
        }
    };
    self.stillDefault = function() {
        return _.isMatch($scope.config, self.default);
    }

    // Project Init
    self.startProject = function() {
        var formData = configFUNC.postForm('check')
        $http({
            method: 'POST', url: configFUNC.connection() + '/query',
            headers: { 'Content-Type': undefined }, data: formData,
            transformRequest: angular.identity
        })
        .success(function (data) {
            console.log(data)
            if (data.exists) {
                swal({ title: "Job Exists",
                       text: "Continuing can overwrite the existing job!",
                       type: "warning", showCancelButton: true,
                       confirmButtonText: "Unlock Job", cancelButtonText: "Change Job Name",
                       closeOnConfirm: false, closeOnCancel: true},
                       function(isConfirm){
                        if (isConfirm) {
                            swal({ title: "Project PIN", text: "Unlock the Project", type: "input",
                                   showCancelButton: true, closeOnConfirm: false, animation: "slide-from-top",
                                   inputPlaceholder: "****"
                                 },
                                 function(inputValue){
                                    if (inputValue === false) return false;
                                    if (inputValue === "") {
                                        swal.showInputError("You need to add a 4 number PIN code.");
                                        return false
                                    }
                                    else {
                                        if (parseInt(inputValue) === parseInt(data.pin.join(''))) {
                                            $http({
                                                method: 'POST', url: configFUNC.connection() + '/query',
                                                headers: { 'Content-Type': undefined }, data: configFUNC.postForm('load'),
                                                transformRequest: angular.identity
                                            })
                                            .success(function (data) {
                                                if (data.success == 'ok') {
                                                    console.log(data)
                                                    self.JSON2project(data.data);
                                                    swal("Nice!", "You unlocked the Project", "success");
                                                }
                                                else {
                                                    swal("Error!", "Unable to load old job", "error");
                                                }
                                            })
                                            .error(function (data, status) {
                                                console.log(data); console.log(status);
                                                swal("Error!", "http error: " + status, "error");
                                            });
                                        }
                                        else {
                                            swal("Error!", "That's not the correct PIN.", "error");
                                            $scope.config.name = null; $scope.$apply();
                                        }
                                    }
                                 }
                                );
                        } else { $scope.config.name = null; $scope.$apply(); }
                       }
                    );
            }
            else if (!data.exists) {
                var s = "000" + Math.floor(Math.random() * 1000)
                s = s.substr(s.length-4)
                $scope.config.pin = _.map(s.substr(s.length-4).split(''), function(num){return parseInt(num)});
                var d = new Date()
                $scope.config.date = [d.getFullYear(), d.getMonth() + 1, d.getDate()];
            }
            else { swal("Error!", "Unknown error has occurred.", "error"); }

        })
        .error(function (data, status) {
            console.log(data); console.log(status);
            swal("Error!", "http error: " + status, "error");
        });
        self.name2default();
    }
});