/*
* @Author: jaume.bonet
* @Date:   2016-11-02 13:27:34
* @Last Modified by:   jaume.bonet
* @Last Modified time: 2016-11-09 18:33:00
*/

'use strict';

var BigEditorCanvas   = new fabric.Canvas('gridcanvas', { selection: false });

fabric.Canvas.prototype.makeGridCanvas = function( grid, do_center, do_axis ) {

    var lineprop = { evented: false, selectable: false, hasBorders: false, hasControls: false };
    // Vertical Rule
    for (var i = 1; i < (this.getWidth() / grid); i++) {
        var line = new fabric.Line([ i * grid, 0, i * grid, this.getHeight()], lineprop);
        this.add(line);
        if ((i * grid % (grid * 5)) !== 0) line.setStroke('#d3d3d3').sendToBack();
        else                               line.setStroke('#3f3f3f').bringToFront();
    }
    // Horizontal Rule
    for (var i = 1; i < (this.getHeight() / grid); i++) {
        var line = new fabric.Line([ 0, i * grid, this.getWidth(), i * grid], lineprop);
        this.add(line);
        if ((i * grid % (grid * 5)) !== 0) line.setStroke('#d3d3d3').sendToBack();
        else                               line.setStroke('#3f3f3f').bringToFront();
    }
    // Center
    if (do_center) {
        var circprop = { radius: grid * 5, evented: false, selectable: false, hasBorders: false, hasControls: false, fill:undefined, stroke:'#e4e4e4', strokeWidth: 1};
        var centre_circle = new fabric.Circle(circprop);
        this.add(centre_circle);
        centre_circle.sendToBack().center();
    }
    // Axis
    if (do_axis) {
        var point = grid * 2.5;
        var line = new fabric.Line([ point, this.getHeight() - point, point, this.getHeight() - point - (grid * 5)], lineprop);
        this.add(line);
        line.setStroke('#000000').setStrokeWidth(2).bringToFront();
        line = new fabric.Line([ point, this.getHeight() - point, point + (grid * 5), this.getHeight() - point], lineprop);
        this.add(line);
        line.setStroke('#000000').setStrokeWidth(2).bringToFront();
        var base   = grid * 1.3;
        var height = (base / 2) * Math.sqrt(3);
        var arrow  = new fabric.Triangle({width: base, height: height, left: point - (base / 2), top: this.getHeight() - point - (grid * 5) - height, evented: false, selectable: false, hasBorders: false, hasControls: false});
        this.add(arrow);
        arrow.setStroke('#000000').setStrokeWidth(2).setFill('#000000').bringToFront();
        arrow  = new fabric.Triangle({width: base, height: height, left: point + (grid * 5) - 1, top: this.getHeight() - point - height + (base / 2) - 0.5, evented: false, selectable: false, hasBorders: false, hasControls: false});
        this.add(arrow);
        arrow.setStroke('#000000').setStrokeWidth(2).setFill('#000000').setAngle(90).bringToFront();
        var text = new fabric.IText('z', { left: point + grid*0.9, top: this.getHeight() - point - (grid * 7.2), fontSize: grid * 2.5, useNative: true, evented: false, selectable: false, hasBorders: false, hasControls: false });
        this.add(text);
        text = new fabric.IText('x', { left: point + (grid * 5.2), top: this.getHeight() - point - grid*3, fontSize: grid * 2.5, useNative: true, evented: false, selectable: false, hasBorders: false, hasControls: false });
        this.add(text);
    }
};

fabric.Canvas.prototype.getEventedItemsStack = function() {
  var objectList = [],
      objects = this.getObjects();

      for (var i = this.size() - 1; i >= 0; i--) {
          if (objects[i].get('evented')) {
            objectList.push(objects[i]);
          }
          else {
            break;
          }
      }
  return objectList;
};

fabric.Canvas.prototype.centerEventedItems = function() {
    var coord = {x: 0, y: 0, c: 0};
    var sse = this.getEventedItemsStack();
    if (sse.length > 0) {
        for (var i = sse.length - 1; i >= 0; i--) {
            var center = sse[i].getCenterPoint();
            coord.x = coord.x + center.x;
            coord.y = coord.y + center.y;
            coord.c = coord.c + 1;
        }
        coord.x = coord.x / coord.c;
        coord.y = coord.y / coord.c;
        for (var i = sse.length - 1; i >= 0; i--) {
            sse[i].setLeft(sse[i].getLeft() - coord.x + this.getCenter().left);
            sse[i].setTop(sse[i].getTop() - coord.y + this.getCenter().top);
            sse[i].setCoords();
        }
        this.renderAll();
    }
};

fabric.Canvas.prototype.lockCanvas = function() {
    var sse = this.getEventedItemsStack();
    for (var i = sse.length - 1; i >= 0; i--) {
        sse[i].set('lockMovementX', true);
        sse[i].set('lockMovementY', true);
        sse[i].set('lockRotation', true);
    }
    return true;
};

fabric.Canvas.prototype.unlockCanvas = function() {
    var sse = this.getEventedItemsStack();
    for (var i = sse.length - 1; i >= 0; i--) {
        sse[i].set('lockMovementX', false);
        sse[i].set('lockMovementY', false);
        sse[i].set('lockRotation', false);
    }
    return false;
};

fabric.Object.prototype.initTopoValues = function(type) {

    this.tpb           = {};
    this.tpb.id        = null;
    this.tpb.type      = type;
    this.tpb.length    = 0;
    this.tpb.ref       = null
    this.tpb.shift     = {x: 0, y: 0, z: 0};
    this.tpb.tilt      = {x: 0, y: 0, z: 0};
    this.tpb.sequence  = {seq: "", editable: true};
    this.tpb.structure = {seq: "", editable: true};
    this.tpb.edge      = 0;
    return this;
}

fabric.Object.prototype.updateTopoValues = function( center , proportion ) {

    this.tpb.tilt.y  = this.angle;
    this.tpb.shift.x = (this.left - center.left) / proportion;
    this.tpb.shift.z = (this.top - center.top) / proportion;

    return this;
}

fabric.Object.prototype.updateValuesTopo = function( center , proportion ) {

    this.angle = this.tpb.tilt.y;
    this.left = (this.tpb.shift.x * proportion) + center.left;
    this.top = (this.tpb.shift.z * proportion) + center.top ;

    return this;
}