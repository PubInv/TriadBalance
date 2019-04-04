// Copyright 2019, Robert L. Read
// This file is part of TriadBalance.
//
// TriadBalance is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TriadBalance is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TriadBalance.  If not, see <https://www.gnu.org/licenses/>.

"use strict";


// This test file is intended to be run in the browser.
// Sadly, I can't figure out how to make exactly the same
// code run from the command line and browser.
// Since this is a GUI element, I am opting for the browser.
// To use mocha from the cli, uncomment these includes.

// var expect = require('chai').expect;
// var should = require('chai').should();

describe('sum', function () {
  it('should return sum of arguments', function () {
    chai.expect(1+2).to.equal(3);
  });
});
