/*
* @Author: jaume.bonet
* @Date:   2016-10-31 17:50:44
* @Last Modified by:   jaume.bonet
* @Last Modified time: 2016-11-10 10:17:59
*/

'use strict';

function makeHelix(x, y, proportion, placeholder) {
    var radius = 2.3 * proportion;
    var stroke = 2;
    var radiusPlus = radius + (stroke / 2);
    var fill = !placeholder ? '#ccccff' : '#fdfefe';  
    var strokec = !placeholder ? '#000000' : '#e6b0aa';
    var circle = new fabric.Circle({ radius: radius, fill:fill, stroke:strokec, strokeWidth: stroke});
    var line   = new fabric.Line([radiusPlus, radiusPlus, radiusPlus, 0], {stroke:strokec, strokeWidth: stroke / 2});
    return new fabric.Group([circle, line], {
        top: y, left: x,
        lockScalingX: true, lockScalingY: true,
        originX: 'center', originY: 'center'}
    ).initTopoValues('H');
}

function makeBeta(x, y, proportion, placeholder) {
    var base   = 20;
    var height = (base / 2) * Math.sqrt(3);
    var stroke = 2;
    var fill = !placeholder ? '#ffcccc' : '#fdfefe';  
    var strokec = !placeholder ? '#000000' : '#e6b0aa';
    return new fabric.Triangle( {
        width: base, height: height, left: x, top: y,
        fill: fill, stroke:strokec, strokeWidth: stroke,
        lockScalingX: true, lockScalingY: true,
        originX: 'center', originY: 'center'}
    ).initTopoValues('E');
}

function makeCross(x, y, proportion, placeholder) {
    var stroke = 2;
    var size   = 20 / 3;
    var fill = !placeholder ? '#ccffcc' : '#fdfefe';  
    var strokec = !placeholder ? '#000000' : '#e6b0aa';

    var crosspoints = [ {x: 0, y: size}, {x: 0, y: size * 2}, {x: size, y: size * 2},
                        {x: size, y: size * 3}, {x: size * 2, y: size * 3},
                        {x: size * 2, y: size * 2}, {x: size * 3, y: size * 2},
                        {x: size * 3, y: size}, {x: size * 2, y: size},
                        {x: size * 2, y: 0}, {x: size, y: 0}, {x: size, y: size}]

    return new fabric.Polygon(crosspoints, {
        left: x, top: y, angle: 45,
        fill: fill, stroke:strokec, strokeWidth: stroke,
        lockScalingX: true, lockScalingY: true,
        originX: 'center', originY: 'center'
      }
    ).initTopoValues('X');
}

function sseMaker(type, center, config) {
    var sse;
    switch(type) {
        case 'helix':
        case 'H':
            sse = makeHelix(center.left, center.top, config.interface.proportion, false);
            sse.tpb.length = config.hLength;
            break;
        case 'beta':
        case 'E':
            sse = makeBeta(center.left, center.top, false);
            sse.tpb.length = config.eLength;
            break;
        case 'cross':
        case 'X':
            sse = makeCross(center.left, center.top, false);
            sse.tpb.length = config.eLength;
            break;
    }
    return sse;
}

function placeholder2sse( placeholder, config ){
    var type = placeholder.tpb.type;
    var center = {left: placeholder.left, top: placeholder.top};
    return sseMaker(type, center, config);
}

function ssePlaceHolders(guide, request, config) {
    var plchldrs = [];
    var dist;
    switch(request) {
        case 'helix':
            if (guide.tpb.type == 'H') dist = config.hhDist * config.interface.proportion;
            if (guide.tpb.type == 'E') dist = config.heDist * config.interface.proportion;
            plchldrs.push(makeHelix(guide.left + dist, guide.top, config.interface.proportion, true))
            plchldrs.push(makeHelix(guide.left - dist, guide.top, config.interface.proportion, true))
            plchldrs.push(makeHelix(guide.left, guide.top + dist, config.interface.proportion, true))
            plchldrs.push(makeHelix(guide.left, guide.top - dist, config.interface.proportion, true))
            break;
        case 'beta':
            if (guide.tpb.type == 'H') dist = config.heDist * config.interface.proportion;
            if (guide.tpb.type == 'E') dist = config.epDist * config.interface.proportion;
            plchldrs.push(makeBeta(guide.left + dist, guide.top, config.interface.proportion, true))
            plchldrs.push(makeBeta(guide.left - dist, guide.top, config.interface.proportion, true))
            if (guide.tpb.type == 'E') dist = config.etDist * config.interface.proportion;
            plchldrs.push(makeBeta(guide.left, guide.top + dist, config.interface.proportion, true))
            plchldrs.push(makeBeta(guide.left, guide.top - dist, config.interface.proportion, true))
            break;
    }
    return plchldrs
}

function makeMotifGroup( mtf, config ){
    var eye = makeBeta( 0, 10 * config.interface.proportion, config.interface.proportion, true);
    eye.setStroke('#000000');
    eye.setFill('#ffffff');
    var g = new fabric.Group([], { lockScalingX: true, lockScalingY: true,
                                   originX: 'center', originY: 'center'});
    for (var i = mtf.length - 1; i >= 0; i--) {
        var left = mtf[i].center[0] * config.interface.proportion;
        var top  = mtf[i].center[2] * config.interface.proportion;
        var str = sseMaker(mtf[i].type, { left: left, top: top}, config);
        g.addWithUpdate(str);
    }
    var h = g.getHeight();
    var w = g.getWidth();
    var r = new fabric.Rect({ left: 0, top: 0, fill:'#ffffff', stroke: '#000000', width: w, height: h, strokeWidth: 1, originX: 'center', originY: 'center' });
    g.add(r);
    r.sendToBack();
    g.add(eye);
    var l = new fabric.Line([0, 0, 0, 10 * config.interface.proportion], {stroke:'#000000', strokeWidth: 2})
    g.add(l);
    l.sendToBack();
    return g;
};