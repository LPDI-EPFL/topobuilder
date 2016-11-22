/*
* @Author: jaume.bonet
* @Date:   2016-11-10 13:19:04
* @Last Modified by:   jaume.bonet
* @Last Modified time: 2016-11-11 14:52:20
*/

'use strict';

angular.module('topobuilder')
.directive('jbConfigManager', function(){
    return {
        scope: false,
        templateUrl:'app/self/static/configmanager.html'
    };
})
.directive('jbCodingManager', function(){
    return {
        scope: false,
        templateUrl:'app/self/static/codingmanager.html'
    };
})
.directive('lowercased', function() {
    return {
        require: 'ngModel',
        link: function(scope, element, attrs, modelCtrl) {
            modelCtrl.$parsers.push(function(input) {
                return input ? input.toLowerCase() : "";
            });
            element.css("text-transform", "lowercase");
        }
    };
})
.directive('reactOnEnter', function () {
    return function (scope, element, attrs) {
        element.bind("keydown keypress", function (event) {
           var key = typeof event.which === "undefined" ? event.keyCode : event.which;
            if(key === 13) {
                var content = scope.$eval(attrs.ngModel);
                if (typeof content != 'undefined' && content) {
                    scope.$apply(function (){
                        scope.$eval(attrs.reactOnEnter);
                    });
                    event.preventDefault();
                }
            }
        });
    };
}).directive('folderSecureText', function() {
    return {
        restrict: 'A',
        link: function($scope, $element) {
            $element.bind('keydown', function(event) {
                var key = typeof event.which === "undefined" ? event.keyCode : event.which;
                if(key != 13 && key != 8) { // Skip Enter key for compatibility with reactOnEnter; also with delete
                    if (event.altKey == true || event.ctrlKey == true) { event.preventDefault(); }
                    if (!((event.keyCode >= 48 && event.keyCode <= 57 && event.shiftKey== false) ||
                          (event.keyCode >= 65 && event.keyCode <= 90) || (event.key == '_') ||
                          (event.keyCode >= 97 && event.keyCode <= 122))) {
                        event.preventDefault();
                        console.log(event)
                    }
                }
            });
        }
    }
}).directive('hostSecureText', function() {
    return {
        restrict: 'A',
        link: function($scope, $element) {
            $element.bind('keydown', function(event) {
                var key = typeof event.which === "undefined" ? event.keyCode : event.which;
                if(key != 13 && key != 8) { // Skip Enter key for compatibility with reactOnEnter; also with delete
                    if (event.altKey == true || event.ctrlKey == true) { event.preventDefault(); }
                    if (!((event.keyCode >= 48 && event.keyCode <= 57 && event.shiftKey== false) ||
                          (event.keyCode >= 65 && event.keyCode <= 90) || (event.key == '.') ||
                          (event.keyCode >= 97 && event.keyCode <= 122))) {
                        event.preventDefault();
                        console.log(event)
                    }
                }
            });
        }
    }
}).directive('focusMe', function($timeout) {
  return {
    link: function(scope, element, attrs) {
        scope.$watch(attrs.focusMe, function(value) {
            if(value === true) { 
              $timeout(function() {
                    element[0].focus();
                    scope[attrs.focusMe] = false;
              });
            }
        });
    }
  };
}).directive('bindFile', function () {
    return {
        require: "ngModel",
        restrict: 'A',
        link: function ($scope, el, attrs, ngModel) {
            el.bind('change', function (event) {
                ngModel.$setViewValue(event.target.files[0]);
                $scope.$apply();
            });
            $scope.$watch(function () {
                return ngModel.$viewValue;
            }, function (value) {
                if (!value) {
                    el.val("");
                }
            });
        }
    };
});